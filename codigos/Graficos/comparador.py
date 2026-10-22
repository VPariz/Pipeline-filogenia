# -*- coding: utf-8 -*-
"""
comparador.py
--------------------
Compara duas arvores .nwk. Calcula a diferenca topologica (RF) e 
nossa metrica customizada (CTS / Diff Prof).

NOTA: O CTS e o p-value calculados aqui sao os usados na versao final do artigo!
"""

from __future__ import annotations

import os
import random

import numpy as np
from scipy import stats

import visualizador_arvores as va
from ete3 import Tree

N_BASELINE_PADRAO = 30


def _carregar_ete3(caminho_nwk: str) -> Tree:
    with open(caminho_nwk, "r", encoding="utf-8", errors="replace") as f:
        texto = f.read().strip()
    t = Tree(texto, format=1)
    for leaf in t.get_leaves():
        if leaf.name.endswith(".txt"):
            leaf.name = leaf.name[:-4]
    return t


def _nos_internos_validos(tree: Tree):
    n = len(tree.get_leaves())
    resultado = []
    for node in tree.traverse():
        if node.is_leaf() or node.is_root():
            continue
        fs = frozenset(l.name for l in node.get_leaves())
        if 1 < len(fs) < n:
            resultado.append((node, fs))
    return resultado


def _canonizar(fs: frozenset, taxa_totais: frozenset) -> frozenset:
    complemento = taxa_totais - fs
    candidatos = sorted([fs, complemento], key=lambda s: (len(s), sorted(s)))
    return candidatos[0]


def analisar_biparticoes(caminho_a: str, caminho_b: str) -> dict:
    ta = _carregar_ete3(caminho_a)
    tb = _carregar_ete3(caminho_b)

    folhas_a = frozenset(l.name for l in ta.get_leaves())
    folhas_b = frozenset(l.name for l in tb.get_leaves())
    taxa_totais = folhas_a | folhas_b

    nos_a = _nos_internos_validos(ta)
    nos_b = _nos_internos_validos(tb)

    bip_a_canon = {_canonizar(fs, taxa_totais) for _, fs in nos_a}
    bip_b_canon = {_canonizar(fs, taxa_totais) for _, fs in nos_b}
    biparticoes_comuns = bip_a_canon & bip_b_canon

    clados_a_raw = {fs for _, fs in nos_a}

    detalhes_b = []
    for _, fs in nos_b:
        detalhes_b.append({
            "frozenset": fs,
            "correto_rf": _canonizar(fs, taxa_totais) in biparticoes_comuns,
            "correto_cts": fs in clados_a_raw,           # antigo "correto_diffprof"
        })

    total_b = len(nos_b)
    corretos_rf = sum(1 for d in detalhes_b if d["correto_rf"])
    corretos_cts = sum(1 for d in detalhes_b if d["correto_cts"])

    rf_absoluto, rf_maximo, nrf, erro_rf = None, None, None, None
    try:
        r = ta.compare(tb, unrooted=True)
        rf_absoluto = r.get("rf")
        rf_maximo = r.get("max_rf")
        nrf = (rf_absoluto / rf_maximo) if (rf_absoluto is not None and rf_maximo) else None
    except Exception as exc:
        erro_rf = str(exc)

    return {
        "arquivo_a": os.path.basename(caminho_a),
        "arquivo_b": os.path.basename(caminho_b),
        "taxa_totais": taxa_totais,
        "biparticoes_comuns": biparticoes_comuns,
        "clados_a_raw": clados_a_raw,
        "detalhes_b": detalhes_b,
        "total_b": total_b,
        "corretos_rf": corretos_rf,
        "corretos_cts": corretos_cts,
        "rf_absoluto": rf_absoluto,
        "rf_maximo": rf_maximo,
        "nrf": nrf,
        "erro_rf": erro_rf,
        "taxons_exclusivos_a": sorted(folhas_a - folhas_b),
        "taxons_exclusivos_b": sorted(folhas_b - folhas_a),
    }


# =========================================================================== #
# Motor de calculo dinamico (v3/v5)
# =========================================================================== #
def alinhar_arvores(t1: Tree, t2: Tree, ocultar_outgroup: bool = True):
    t1c = t1.copy()
    t2c = t2.copy()

    for t in (t1c, t2c):
        try:
            t.resolve_polytomy(recursive=True)
        except Exception:
            pass

    folhas1 = {l.name for l in t1c.get_leaves()}
    folhas2 = {l.name for l in t2c.get_leaves()}
    comuns = folhas1 & folhas2

    if ocultar_outgroup:
        comuns = {n for n in comuns if not va.eh_outgroup_artificial(n)}

    comuns = sorted(comuns)
    if len(comuns) < 3:
        raise ValueError(
            f"Apenas {len(comuns)} táxon(ns) em comum (mínimo: 3)."
        )

    t1c.prune(comuns, preserve_branch_length=True)
    t2c.prune(comuns, preserve_branch_length=True)
    return t1c, t2c, comuns


def get_patristic_vector(tree: Tree, sorted_leaves: list) -> np.ndarray:
    n = len(sorted_leaves)
    valores = np.empty(n * (n - 1) // 2, dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            valores[idx] = tree.get_distance(sorted_leaves[i], sorted_leaves[j])
            idx += 1
    return valores


def calcular_fit_distancias(vec_ref: np.ndarray, vec_target: np.ndarray) -> dict:
    if len(vec_ref) < 2 or len(vec_target) < 2 or len(vec_ref) != len(vec_target):
        return {"fator_escala": None, "pearson_r": None, "spearman_r": None, "erro_escala": None}

    resultado = {"fator_escala": None, "pearson_r": None, "spearman_r": None, "erro_escala": None}

    denom = float(np.dot(vec_target, vec_target))
    if denom > 0:
        fator_escala = float(np.dot(vec_ref, vec_target) / denom)
        resultado["fator_escala"] = fator_escala
        norma_ref = float(np.linalg.norm(vec_ref))
        if norma_ref > 0:
            residuo = vec_ref - fator_escala * vec_target
            resultado["erro_escala"] = float(np.linalg.norm(residuo) / norma_ref)

    if np.std(vec_ref) > 0 and np.std(vec_target) > 0:
        try:
            resultado["pearson_r"] = float(stats.pearsonr(vec_ref, vec_target)[0])
        except Exception:
            pass
        try:
            resultado["spearman_r"] = float(stats.spearmanr(vec_ref, vec_target)[0])
        except Exception:
            pass

    return resultado


def _alturas_topologicas(tree: Tree) -> dict:
    alturas = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            alturas[node] = 0
        else:
            alturas[node] = 1 + max(alturas[c] for c in node.children)
    return alturas


def calcular_cts(t_ref_ete: Tree, t_tar_ete: Tree) -> float:
    """
    Calcula o Custom Topological Score (CTS) - antigo Diff Prof.
    Basicamente, olha pros clados e da um peso maior pra erros na raiz da arvore.
    """
    alturas_ref = _alturas_topologicas(t_ref_ete)
    alturas_tar = _alturas_topologicas(t_tar_ete)

    clados_ref = {
        frozenset(l.name for l in n.get_leaves()): alturas_ref[n]
        for n in t_ref_ete.traverse()
        if not n.is_leaf() and not n.is_root()
    }

    total = 0.0
    contagem = 0
    for n in t_tar_ete.traverse():
        if n.is_leaf() or n.is_root():
            continue
        contagem += 1
        fs = frozenset(l.name for l in n.get_leaves())
        if fs in clados_ref:
            total += abs(alturas_tar[n] - clados_ref[fs])
        else:
            total += alturas_tar[n]

    return (total / contagem) if contagem else 0.0


def _pvalue_topologia(t_ref: Tree, t_tar: Tree, nrf_observado, n_baseline: int):
    if nrf_observado is None:
        return None

    nomes_originais = [l.name for l in t_tar.get_leaves()]
    contagem, validos = 0, 0

    for _ in range(n_baseline):
        tar_copia = t_tar.copy()
        folhas_copia = tar_copia.get_leaves()
        permutados = nomes_originais[:]
        random.shuffle(permutados)
        for folha, novo_nome in zip(folhas_copia, permutados):
            folha.name = novo_nome
        try:
            r = t_ref.compare(tar_copia, unrooted=True)
            rf, max_rf = r.get("rf"), r.get("max_rf")
            if max_rf:
                validos += 1
                if (rf / max_rf) <= nrf_observado:
                    contagem += 1
        except Exception:
            continue

    return (contagem / validos) if validos else None


def _pvalue_distancias(vec_ref: np.ndarray, vec_target: np.ndarray,
                        pearson_observado, n_baseline: int):
    if pearson_observado is None or len(vec_target) < 2:
        return None

    contagem, validos = 0, 0
    indices = np.arange(len(vec_target))
    for _ in range(n_baseline):
        permutado = vec_target[np.random.permutation(indices)]
        if np.std(permutado) == 0:
            continue
        try:
            r_perm = stats.pearsonr(vec_ref, permutado)[0]
        except Exception:
            continue
        validos += 1
        if abs(r_perm) >= abs(pearson_observado):
            contagem += 1

    return (contagem / validos) if validos else None


def executar_comparacao_completa(caminho_ref: str, caminho_tar: str,
                                  calcular_pvalue: bool = False,
                                  n_baseline: int = N_BASELINE_PADRAO) -> dict:
    resultado = {
        "arquivo_ref": os.path.basename(caminho_ref),
        "arquivo_tar": os.path.basename(caminho_tar),
        "taxons": None,
        "rf_absoluto": None, "rf_maximo": None, "nrf": None,
        "cts": None,                      # antigo "diff_prof"
        "pearson_r": None, "spearman_r": None, "erro_escala": None,
        "p_value_topologia": None, "p_value_distancias": None,
        "erro": None,
    }

    try:
        t_ref = _carregar_ete3(caminho_ref)
        t_tar = _carregar_ete3(caminho_tar)
        t_ref_a, t_tar_a, taxa_comuns = alinhar_arvores(t_ref, t_tar)
        resultado["taxons"] = len(taxa_comuns)

        r = t_ref_a.compare(t_tar_a, unrooted=True)
        rf_absoluto, rf_maximo = r.get("rf"), r.get("max_rf")
        resultado["rf_absoluto"] = rf_absoluto
        resultado["rf_maximo"] = rf_maximo
        resultado["nrf"] = (rf_absoluto / rf_maximo) if (rf_absoluto is not None and rf_maximo) else None

        resultado["cts"] = calcular_cts(t_ref_a, t_tar_a)   # funcao renomeada

        vec_ref = get_patristic_vector(t_ref_a, taxa_comuns)
        vec_tar = get_patristic_vector(t_tar_a, taxa_comuns)
        fit = calcular_fit_distancias(vec_ref, vec_tar)
        resultado["pearson_r"] = fit["pearson_r"]
        resultado["spearman_r"] = fit["spearman_r"]
        resultado["erro_escala"] = fit["erro_escala"]

        if calcular_pvalue:
            resultado["p_value_topologia"] = _pvalue_topologia(
                t_ref_a, t_tar_a, resultado["nrf"], n_baseline)
            resultado["p_value_distancias"] = _pvalue_distancias(
                vec_ref, vec_tar, resultado["pearson_r"], n_baseline)

    except Exception as exc:
        resultado["erro"] = str(exc)

    return resultado


# --------------------------------------------------------------------------- #
# Geradores de color_fn (v5: nomes atualizados)
# --------------------------------------------------------------------------- #
def gerador_cor_rf(analise: dict):
    comuns = analise["biparticoes_comuns"]
    taxa_totais = analise["taxa_totais"]
    n_total = len(taxa_totais)

    def color_fn(fs: frozenset):
        if len(fs) <= 1 or len(fs) >= n_total:
            return {"color": va.COR_RAMO_PADRAO, "linestyle": "-", "linewidth": 1.5}
        fs_canon = _canonizar(fs, taxa_totais)
        if fs_canon in comuns:
            return {"color": "#2ecc71", "linestyle": "-", "linewidth": 2.5}
        return {"color": "#e74c3c", "linestyle": "--", "linewidth": 2.0}

    return color_fn


def gerador_cor_cts(analise: dict):   # antigo gerador_cor_diffprof
    clados_a = analise["clados_a_raw"]
    n_total = len(analise["taxa_totais"])

    def color_fn(fs: frozenset):
        if len(fs) <= 1 or len(fs) >= n_total:
            return {"color": va.COR_RAMO_PADRAO, "linestyle": "-", "linewidth": 1.5}
        if fs in clados_a:
            return {"color": "#27ae60", "linestyle": "-", "linewidth": 2.5}
        return {"color": "#c0392b", "linestyle": "--", "linewidth": 2.0}

    return color_fn


# --------------------------------------------------------------------------- #
# Compatibilidade com v1
# --------------------------------------------------------------------------- #
def comparar_arvores(caminho_a: str, caminho_b: str) -> dict:
    analise = analisar_biparticoes(caminho_a, caminho_b)
    return {
        "arquivo_a": analise["arquivo_a"],
        "arquivo_b": analise["arquivo_b"],
        "taxons_exclusivos_a": analise["taxons_exclusivos_a"],
        "taxons_exclusivos_b": analise["taxons_exclusivos_b"],
        "rf_absoluto": analise["rf_absoluto"],
        "rf_maximo": analise["rf_maximo"],
        "nrf": analise["nrf"],
        "erro": analise["erro_rf"],
    }


if __name__ == "__main__":
    print("comparador.py v5 carregado com sucesso (CTS).")
# -*- coding: utf-8 -*-
"""
visualizador_arvores.py
------------------------------
Script de plotagem com Matplotlib + Bio.Phylo.
Ele limpa os nomes dos taxons e ajusta a topologia pra renderizar 
arvores ou cladogramas no formato de publicacao do artigo
"""

from __future__ import annotations

import io
import os
import re
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_CGI_STUB_PATH = os.path.join(_THIS_DIR, "cgi.py")
if not os.path.exists(_CGI_STUB_PATH):
    with open(_CGI_STUB_PATH, "w", encoding="utf-8") as _f:
        _f.write('"""Polyfill: cgi removido no Python 3.13"""\n')
sys.path.insert(0, _THIS_DIR)

from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTGROUP_ARTIFICIAL = "randomFasta-N16k"
_INNER_NAME_RE = re.compile(r"^Inner\d+$", re.IGNORECASE)

COR_RAMO_PADRAO = "#333333"
COR_OUTGROUP = "darkorange"


def _nome_limpo(nome: str | None) -> str | None:
    if nome is None:
        return None
    if _INNER_NAME_RE.match(nome):
        return None
    if nome.endswith(".txt"):
        return nome[:-4]
    return nome


def carregar_arvore(caminho_nwk: str):
    with open(caminho_nwk, "r", encoding="utf-8", errors="replace") as f:
        texto = f.read().strip()
    tree = Phylo.read(io.StringIO(texto), "newick")
    for clade in tree.find_clades():
        clade.name = _nome_limpo(clade.name)
    return tree


def eh_outgroup_artificial(nome_folha: str) -> bool:
    return nome_folha is not None and OUTGROUP_ARTIFICIAL in nome_folha


def _prunar_outgroup(tree):
    alvo = None
    for leaf in tree.get_terminals():
        if eh_outgroup_artificial(leaf.name):
            alvo = leaf
            break
    if alvo is not None:
        tree.prune(alvo)


def _converter_para_cladograma(tree):
    """
    Ignora as distancias brutas da metrica NCD e ajusta a estrutura do cladograma
    pra alinhar as folhas a direita e garantir que os ramos fiquem legiveis
    no formato final do artigo
    """
    def calc_height(clade):
        if clade.is_terminal():
            clade._height = 0
            return 0
        h = 1 + max(calc_height(c) for c in clade.clades)
        clade._height = h
        return h

    calc_height(tree.root)
    tree.root.branch_length = 0

    def set_lengths(clade):
        for child in clade.clades:
            child.branch_length = clade._height - child._height
            set_lengths(child)

    set_lengths(tree.root)


def desenhar_cladograma(
    tree,
    ax,
    titulo: str = "",
    ocultar_outgroup: bool = False,
    usar_distancias: bool = True,
    modo_publicacao: bool = False,
    color_fn=None,
    titulo_extra: str = "",
):
    import copy
    t = copy.deepcopy(tree)

    if ocultar_outgroup:
        _prunar_outgroup(t)

    if not usar_distancias:
        _converter_para_cladograma(t)

    if modo_publicacao:
        plt.rcParams["font.family"] = "DejaVu Serif"
        fontSize = 10
    else:
        plt.rcParams["font.family"] = "sans-serif"
        fontSize = 9

    def label_func(clade):
        if clade.is_terminal():
            return f"  {clade.name}"
        return ""

    if color_fn is not None:
        for clade in t.find_clades():
            leaves = frozenset(c.name for c in clade.get_terminals() if c.name)
            estilo = color_fn(leaves)
            clade.color = estilo.get("color", COR_RAMO_PADRAO)

    Phylo.draw(
        t,
        axes=ax,
        do_show=False,
        label_func=label_func,
        branch_labels=None,
    )

    titulo_final = titulo
    if titulo_extra:
        titulo_final = f"{titulo}\n({titulo_extra})" if titulo else titulo_extra

    if titulo_final:
        ax.set_title(titulo_final, fontsize=fontSize + 2, fontweight="bold", pad=12)


    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(usar_distancias)

    if not usar_distancias:
        ax.get_xaxis().set_visible(False)

    ax.yaxis.set_visible(False)


def exportar_figura(fig, caminho_sem_ext: str):
    fig.savefig(f"{caminho_sem_ext}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{caminho_sem_ext}.svg", format="svg", bbox_inches="tight")


# Alias de compatibilidade
desenhar_arvore = desenhar_cladograma
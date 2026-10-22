# -*- coding: utf-8 -*-
"""
parser_metricas.py
------------------------
Script utilitario pra extrair infos do nome dos arquivos .nwk
e ler o txt do relatorio de metricas gerado em batch
"""

from __future__ import annotations

import os
import re
import pandas as pd

# --------------------------------------------------------------------------- #
# Lista de compressores e algoritmos padroes
# --------------------------------------------------------------------------- #
COMPRESSORES_PADRAO = ["bzip2", "zstd", "brotli", "ppmd", "ficd", "aficd"]

ALGORITMOS_MAP = {
    "upgma": "UPGMA",
    "complete_linkage": "Complete Linkage",
    "neighbor_joining": "Neighbor-Joining",
}

_EXPERIMENTAL_KEYWORDS = ["replica_bat", "replica", "asym", "gzip", "lzma", "zlib"]

# --------------------------------------------------------------------------- #
# Funcao pra extrair meta-infos so olhando pro nome do arquivo NWK
# --------------------------------------------------------------------------- #
def parse_nome_arquivo(nome_arquivo: str) -> dict:
    base = os.path.basename(nome_arquivo)
    stem = re.sub(r"\.nwk$", "", base, flags=re.IGNORECASE)

    if stem.lower() == "modelo":
        return {
            "compressor": "—",
            "abordagem": "Modelo (referência)",
            "algoritmo": "—",
            "sem_distancia": False,
            "legado": False,
        }
    if stem.lower() == "artigo":
        return {
            "compressor": "—",
            "abordagem": "Artigo (referência)",
            "algoritmo": "—",
            "sem_distancia": False,
            "legado": False,
        }

    legado = any(kw in stem.lower() for kw in _EXPERIMENTAL_KEYWORDS)

    compressor = "Desconhecido"
    algoritmo = "Desconhecido"
    abordagem = ""

    for c in COMPRESSORES_PADRAO:
        if c in stem.lower():
            compressor = c
            if c == "ficd":
                compressor = "IRR-C"
            elif c == "aficd":
                compressor = "IRR-C/aFicd"
            break

    for key, label in ALGORITMOS_MAP.items():
        if key in stem.lower():
            algoritmo = label
            break

    if compressor == "Desconhecido" and not legado:
        partes = stem.split("_")
        if partes:
            compressor = partes[0]

    if algoritmo == "Desconhecido" and not legado:
        partes = stem.split("_")
        if len(partes) > 1:
            ultimo = partes[-1]
            if ultimo in ALGORITMOS_MAP:
                algoritmo = ALGORITMOS_MAP[ultimo]

    return {
        "compressor": compressor,
        "abordagem": abordagem,
        "algoritmo": algoritmo,
        "sem_distancia": False,
        "legado": legado,
    }


# --------------------------------------------------------------------------- #
# Leitor do txt de relatorio (pega o Diff Prof antigo ou CTS novo)
# --------------------------------------------------------------------------- #
_BLOCK_RE = re.compile(r"ARQUIVO:\s*(?P<arquivo>\S+)\s*\(T[aá]xons:\s*(?P<taxons>\d+)\)")
_FIELD_PATTERNS = {
    "rf_absoluto": re.compile(r"RF Absoluto:\s*([\-\d.]+)"),
    "nrf": re.compile(r"nRF:\s*([\-\d.]+)"),
    # Novo: captura "Diff Prof" ou "CTS" - chave sera "cts"
    "cts": re.compile(r"(?:Diff\s+Prof|CTS):\s*([\-\d.]+)", re.IGNORECASE),
    "pearson_r": re.compile(r"Pearson r:\s*([\-\d.]+)"),
    "spearman_r": re.compile(r"Spearman r:\s*([\-\d.]+)"),
    "erro_escala": re.compile(r"Erro Escala:\s*([\-\d.]+)"),
}
_PVALUE_RE = re.compile(r"p-value:\s*([\-\d.]+)\s*(\(?\*?Significativo\)?)?", re.IGNORECASE)


def _parse_bloco(texto_bloco: str) -> dict:
    dados = {}
    for campo, regex in _FIELD_PATTERNS.items():
        m = regex.search(texto_bloco)
        dados[campo] = float(m.group(1)) if m else None

    pvalues = _PVALUE_RE.findall(texto_bloco)
    dados["p_value_topologia"] = float(pvalues[0][0]) if len(pvalues) >= 1 else None
    dados["p_value_distancias"] = float(pvalues[1][0]) if len(pvalues) >= 2 else None
    dados["topologia_significativa"] = dados["p_value_topologia"] is not None and dados["p_value_topologia"] < 0.05
    dados["distancias_significativa"] = dados["p_value_distancias"] is not None and dados["p_value_distancias"] < 0.05
    return dados


def parse_relatorio(caminho_relatorio: str) -> pd.DataFrame:
    with open(caminho_relatorio, "r", encoding="utf-8", errors="replace") as f:
        conteudo = f.read()

    matches = list(_BLOCK_RE.finditer(conteudo))
    indices = [m.start() for m in matches]
    indices.append(len(conteudo))

    linhas = []
    for i, m in enumerate(matches):
        bloco = conteudo[m.start(): indices[i + 1]]
        arquivo = m.group("arquivo")
        taxons = int(m.group("taxons"))

        metricas = _parse_bloco(bloco)
        meta = parse_nome_arquivo(arquivo)

        linha = {
            "arquivo": arquivo,
            "taxons": taxons,
            **meta,
            **metricas,
        }
        linhas.append(linha)

    df = pd.DataFrame(linhas)
    if not df.empty:
        df = df.sort_values(
            ["cts", "nrf"], ascending=[True, True], na_position="last"
        ).reset_index(drop=True)
    return df


# --------------------------------------------------------------------------- #
# Geracao dinamica (caso a gente mude a arvore base de referencia)
# --------------------------------------------------------------------------- #
def calcular_relatorio_dinamico(caminho_ref_nwk: str, pasta_arvores: str,
                                 calcular_pvalue: bool = False, n_baseline: int = 30,
                                 incluir_experimentais: bool = False,
                                 progress_callback=None) -> pd.DataFrame:
    import comparador as cmp

    arquivos = []
    for f in os.listdir(pasta_arvores):
        if f.lower().endswith(".nwk"):
            arquivos.append((os.path.join(pasta_arvores, f), False))

    if incluir_experimentais:
        pasta_exp = os.path.join(os.path.dirname(pasta_arvores), "experimentais")
        if os.path.isdir(pasta_exp):
            for f in os.listdir(pasta_exp):
                if f.lower().endswith(".nwk"):
                    arquivos.append((os.path.join(pasta_exp, f), True))

    caminho_ref_abs = os.path.abspath(caminho_ref_nwk)
    total = len(arquivos)

    linhas = []
    for i, (caminho_tar, is_exp) in enumerate(arquivos, start=1):
        if os.path.abspath(caminho_tar) == caminho_ref_abs:
            if progress_callback:
                progress_callback(i, total, os.path.basename(caminho_tar))
            continue

        resultado = cmp.executar_comparacao_completa(
            caminho_ref_nwk, caminho_tar,
            calcular_pvalue=calcular_pvalue, n_baseline=n_baseline,
        )

        nome_arquivo = os.path.basename(caminho_tar)
        meta = parse_nome_arquivo(nome_arquivo)
        if is_exp:
            meta["legado"] = True

        p_topo = resultado.get("p_value_topologia")
        p_dist = resultado.get("p_value_distancias")

        linha = {
            "arquivo": nome_arquivo,
            "taxons": resultado.get("taxons"),
            **meta,
            "rf_absoluto": resultado.get("rf_absoluto"),
            "nrf": resultado.get("nrf"),
            "cts": resultado.get("cts"),           # chave unificada
            "pearson_r": resultado.get("pearson_r"),
            "spearman_r": resultado.get("spearman_r"),
            "erro_escala": resultado.get("erro_escala"),
            "p_value_topologia": p_topo,
            "p_value_distancias": p_dist,
            "topologia_significativa": p_topo is not None and p_topo < 0.05,
            "distancias_significativa": p_dist is not None and p_dist < 0.05,
            "erro_calculo": resultado.get("erro"),
        }
        linhas.append(linha)

        if progress_callback:
            progress_callback(i, total, nome_arquivo)

    df = pd.DataFrame(linhas)
    if not df.empty:
        df = df.sort_values(
            ["cts", "nrf"], ascending=[True, True], na_position="last"
        ).reset_index(drop=True)
    return df


def formatar_relatorio_texto(df: pd.DataFrame, referencia: str = "modelo") -> str:
    def fmt(v):
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return "N/A"
        try:
            return f"{float(v):.4f}"
        except (TypeError, ValueError):
            return str(v)

    linhas = [
        f"# Relatório gerado pela ferramenta (v5) -- referência: {referencia}",
        "# ATENÇÃO: se esta referência não for o modelo taxonômico padrão, os",
        "# valores de CTS e p-value são aproximações calculadas em tempo real.",
        "",
    ]
    for _, row in df.iterrows():
        linhas.append(f"ARQUIVO: {row['arquivo']} (Táxons: {row.get('taxons', '?')})")
        linhas.append("  [Topologia]")
        linhas.append(f"    - RF Absoluto: {fmt(row.get('rf_absoluto'))}")
        linhas.append(f"    - nRF:         {fmt(row.get('nrf'))}")
        linhas.append(f"    - CTS:         {fmt(row.get('cts'))}")      # nome atualizado
        linhas.append(f"    - p-value:     {fmt(row.get('p_value_topologia'))}")
        linhas.append("  [Distâncias/Comprimentos]")
        linhas.append(f"    - Pearson r:   {fmt(row.get('pearson_r'))}")
        linhas.append(f"    - Spearman r:  {fmt(row.get('spearman_r'))}")
        linhas.append(f"    - Erro Escala: {fmt(row.get('erro_escala'))}")
        linhas.append(f"    - p-value:     {fmt(row.get('p_value_distancias'))}")
        linhas.append("")
    return "\n".join(linhas)


# --------------------------------------------------------------------------- #
# Logica de cores e status pras tabelas de resultados
# --------------------------------------------------------------------------- #
def classificar_linha(row) -> str:
    if row.get("legado"):
        return "legado"

    cts = row.get("cts")
    nrf = row.get("nrf")

    if cts is None or nrf is None:
        return "ruim"

    if cts <= 3.0 and nrf <= 0.12:
        return "boa"
    if cts <= 5.0 or nrf <= 0.18:
        return "media"
    return "ruim"


# --------------------------------------------------------------------------- #
# Validacao basica da leitura do relatorio
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    exemplo = """
ARQUIVO: bzip2_neighbor_joining.nwk (Táxons: 20)
  [Topologia]
    - RF Absoluto: 2
    - nRF:         0.1176
    - CTS:         5.0
    - p-value:     0.0100
  [Distâncias/Comprimentos]
    - Pearson r:   0.91
    - Spearman r:  0.92
    - Erro Escala: 0.20
    - p-value:     0.4000

ARQUIVO: aficd_neighbor_joining.nwk (Táxons: 20)
  [Topologia]
    - RF Absoluto: 6
    - nRF:         0.1176
    - Diff Prof:   4.5
    - p-value:     0.0300
  [Distâncias/Comprimentos]
    - Pearson r:   0.80
    - Spearman r:  0.81
    - Erro Escala: 0.30
    - p-value:     0.10
"""
    with open("/tmp/relatorio_v5.txt", "w", encoding="utf-8") as f:
        f.write(exemplo)

    df = parse_relatorio("/tmp/relatorio_v5.txt")
    df["classificacao"] = df.apply(classificar_linha, axis=1)
    print(df[["arquivo", "compressor", "algoritmo", "cts", "nrf", "classificacao"]].to_string())
    print("\nOK: parser v5 com retrocompatibilidade Diff Prof/CTS.")
# -*- coding: utf-8 -*-
"""
gerar_figuras_artigo.py
------------------------
Script pra gerar as imagens em alta definicao pro artigo do BSB 2026.

Gera as seguintes imagens:
  - mcc_tree.png           : Nossa arvore base (VertLife)
  - melhor_zstd.png        : A que deu o melhor resultado (zstd + CL)
  - comparacao_ficd_bzip2.png : Ficd vs bzip2 (pra mostrar q o Ficd acerta na raiz)

Uso:
    python gerar_figuras_artigo.py
    python gerar_figuras_artigo.py --destino "C:/algum/lugar/"
"""

from __future__ import annotations

import argparse
import os
import sys

# Força UTF-8 na saída padrão (necessário no Windows com terminal cp1252)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# --- Garante que o módulo encontra os helpers locais ---
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(_THIS_DIR)
sys.path.insert(0, _THIS_DIR)

# Polyfill para Python 3.13 (remove dependência de cgi)
_CGI_STUB = os.path.join(_THIS_DIR, "cgi.py")
if not os.path.exists(_CGI_STUB):
    with open(_CGI_STUB, "w", encoding="utf-8") as _f:
        _f.write('"""Polyfill: cgi removido no Python 3.13"""\n')

import matplotlib
matplotlib.use("Agg")  # Backend sem janela (geração headless)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import visualizador_arvores as va
import comparador as cmp

# ---------------------------------------------------------------------------
# Caminhos padrão
# ---------------------------------------------------------------------------
_PIBIC_ROOT = os.path.abspath(os.path.join(_THIS_DIR, "..", ".."))

PASTA_ARVORES = os.path.join(_PIBIC_ROOT, "Resultados", "arvores")
DESTINO_PADRAO = os.path.join(
    _PIBIC_ROOT, "ProduçãoEscrita", "BSB2026", "rascunho"
)

ARVORE_MODELO = os.path.join(PASTA_ARVORES, "modelo.nwk")
ARVORE_ZSTD_CL = os.path.join(PASTA_ARVORES, "zstd_complete_linkage.nwk")
ARVORE_FICD_NJ = os.path.join(PASTA_ARVORES, "ficd_neighbor_joining.nwk")
ARVORE_BZIP2_NJ = os.path.join(PASTA_ARVORES, "bzip2_neighbor_joining.nwk")

# ---------------------------------------------------------------------------
# Configurações visuais de publicação (SBC / Overleaf)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

DPI = 300
COR_CORRETO = "#27ae60"   # verde — clado correto (CTS)
COR_ERRADO  = "#c0392b"   # vermelho — clado incorreto (CTS)
COR_PADRAO  = "#333333"   # cinza escuro — ramos neutros


# ---------------------------------------------------------------------------
# Figura 1: Árvore MCC (padrão-ouro)
# ---------------------------------------------------------------------------
def gerar_mcc_tree(destino: str, lang: str = "pt", sufixo: str = ""):
    print(f">> Gerando mcc_tree{sufixo}.png ({lang.upper()}) ...")
    tree = va.carregar_arvore(ARVORE_MODELO)

    fig, ax = plt.subplots(figsize=(7.0, 5.5))

    titulo = (
        "Reference MCC Tree (VertLife · 31 nuclear genes · BEAST 2)"
        if lang == "en" else
        "Árvore MCC de Referência (VertLife · 31 genes nucleares · BEAST 2)"
    )

    va.desenhar_cladograma(
        tree, ax,
        titulo="",
        ocultar_outgroup=True,
        usar_distancias=False,
        modo_publicacao=True,
    )
    ax.set_xlabel("")

    saida = os.path.join(destino, f"mcc_tree{sufixo}")
    va.exportar_figura(fig, saida)
    plt.close(fig)
    print(f"   Salvo: {saida}.png  |  {saida}.svg")


# ---------------------------------------------------------------------------
# Figura 2: Melhor configuração – zstd + Ligação Completa
# ---------------------------------------------------------------------------
def gerar_melhor_zstd(destino: str, lang: str = "pt", sufixo: str = ""):
    print(f">> Gerando melhor_zstd{sufixo}.png ({lang.upper()}) ...")

    tree_modelo = va.carregar_arvore(ARVORE_MODELO)
    tree_zstd   = va.carregar_arvore(ARVORE_ZSTD_CL)

    analise = cmp.analisar_biparticoes(ARVORE_MODELO, ARVORE_ZSTD_CL)
    color_fn = cmp.gerador_cor_cts(analise)

    fig, ax = plt.subplots(figsize=(7.0, 5.5))

    tit_zstd = (
        "zstd + Complete Linkage  (CTS = 3.0 · nRF = 0.1176)"
        if lang == "en" else
        "zstd + Ligação Completa  (CTS = 3,0 · nRF = 0,1176)"
    )
    leg_correto = "Recovered clade (CTS)" if lang == "en" else "Clado recuperado (CTS)"
    leg_errado  = "Divergent clade (CTS)" if lang == "en" else "Clado divergente (CTS)"

    # Plotar apenas zstd+CL com coloração por CTS
    va.desenhar_cladograma(
        tree_zstd, ax,
        titulo=tit_zstd,
        ocultar_outgroup=True,
        usar_distancias=False,
        modo_publicacao=True,
        color_fn=color_fn,
    )

    # Legenda manual
    legenda = [
        mpatches.Patch(color=COR_CORRETO, label=leg_correto),
        mpatches.Patch(color=COR_ERRADO,  label=leg_errado),
    ]
    fig.legend(
        handles=legenda,
        loc="lower center",
        ncol=2,
        fontsize=9,
        frameon=True,
        framealpha=0.9,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.tight_layout()

    saida = os.path.join(destino, f"melhor_zstd{sufixo}")
    va.exportar_figura(fig, saida)
    plt.close(fig)
    print(f"   Salvo: {saida}.png  |  {saida}.svg")


# ---------------------------------------------------------------------------
# Figura 3: Painel comparativo Ficd+NJ vs bzip2+NJ (coloração por CTS)
# ---------------------------------------------------------------------------
def gerar_figuras_irrc_bzip2(destino: str, lang: str = "pt", sufixo: str = ""):
    print(f">> Gerando irrc_nj{sufixo}.png e bzip2_nj{sufixo}.png ({lang.upper()}) ...")

    tree_ficd    = va.carregar_arvore(ARVORE_FICD_NJ)
    tree_bzip2   = va.carregar_arvore(ARVORE_BZIP2_NJ)

    analise_ficd  = cmp.analisar_biparticoes(ARVORE_MODELO, ARVORE_FICD_NJ)
    analise_bzip2 = cmp.analisar_biparticoes(ARVORE_MODELO, ARVORE_BZIP2_NJ)

    color_fn_ficd  = cmp.gerador_cor_cts(analise_ficd)
    color_fn_bzip2 = cmp.gerador_cor_cts(analise_bzip2)

    tit_irrc = (
        "IRR-C + NJ\n(CTS = 4.5 · nRF = 0.2941)"
        if lang == "en" else
        "IRR-C + NJ\n(CTS = 4,5 · nRF = 0,2941)"
    )
    tit_bzip2 = (
        "bzip2 + NJ\n(CTS = 5.0 · nRF = 0.1176)"
        if lang == "en" else
        "bzip2 + NJ\n(CTS = 5,0 · nRF = 0,1176)"
    )
    leg_correto = "Correctly recovered clade (CTS)" if lang == "en" else "Clado recuperado corretamente (CTS)"
    leg_errado  = "Divergent clade from gold standard (CTS)" if lang == "en" else "Clado divergente do padrão-ouro (CTS)"

    legenda = [
        mpatches.Patch(color=COR_CORRETO, label=leg_correto),
        mpatches.Patch(color=COR_ERRADO,  label=leg_errado),
    ]

    # Figura IRR-C + NJ
    fig1, ax1 = plt.subplots(figsize=(7.0, 5.5))
    va.desenhar_cladograma(
        tree_ficd, ax1,
        titulo=tit_irrc,
        ocultar_outgroup=True,
        usar_distancias=False,
        modo_publicacao=True,
        color_fn=color_fn_ficd,
    )
    fig1.legend(handles=legenda, loc="lower center", ncol=2, fontsize=9, frameon=True, framealpha=0.9, bbox_to_anchor=(0.5, -0.02))
    fig1.tight_layout()
    saida1 = os.path.join(destino, f"irrc_nj{sufixo}")
    va.exportar_figura(fig1, saida1)
    plt.close(fig1)

    # Figura bzip2 + NJ
    fig2, ax2 = plt.subplots(figsize=(7.0, 5.5))
    va.desenhar_cladograma(
        tree_bzip2, ax2,
        titulo=tit_bzip2,
        ocultar_outgroup=True,
        usar_distancias=False,
        modo_publicacao=True,
        color_fn=color_fn_bzip2,
    )
    fig2.legend(handles=legenda, loc="lower center", ncol=2, fontsize=9, frameon=True, framealpha=0.9, bbox_to_anchor=(0.5, -0.02))
    fig2.tight_layout()
    saida2 = os.path.join(destino, f"bzip2_nj{sufixo}")
    va.exportar_figura(fig2, saida2)
    plt.close(fig2)

    print(f"   Salvo: {saida1}.png  |  {saida2}.png")


# ---------------------------------------------------------------------------
# Ponto de entrada
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Gera as figuras do artigo BSB 2026 a partir das árvores em Resultados/arvores/"
    )
    parser.add_argument(
        "--destino", "-d",
        default=DESTINO_PADRAO,
        help=f"Pasta de saída das imagens (padrão: {DESTINO_PADRAO})"
    )
    parser.add_argument(
        "--lang", "-l",
        choices=["pt", "en", "both"],
        default="pt",
        help="Idioma das figuras: 'pt' (português), 'en' (inglês) ou 'both' (gera ambos os idiomas)"
    )
    args = parser.parse_args()

    destino = os.path.abspath(args.destino)
    os.makedirs(destino, exist_ok=True)

    print(f"\nPasta de arvores : {PASTA_ARVORES}")
    print(f"Destino das figuras: {destino}")
    print(f"Idioma escolhido : {args.lang.upper()}\n")

    # Verifica se os arquivos necessários existem
    necessarios = {
        "modelo.nwk":                ARVORE_MODELO,
        "zstd_complete_linkage.nwk": ARVORE_ZSTD_CL,
        "ficd_neighbor_joining.nwk": ARVORE_FICD_NJ,
        "bzip2_neighbor_joining.nwk": ARVORE_BZIP2_NJ,
    }
    faltando = [nome for nome, caminho in necessarios.items() if not os.path.exists(caminho)]
    if faltando:
        print("ERRO: Os seguintes arquivos não foram encontrados em Resultados/arvores/:")
        for f in faltando:
            print(f"   - {f}")
        sys.exit(1)

    idiomas = ["pt", "en"] if args.lang == "both" else [args.lang]

    for idioma in idiomas:
        sufixo = "" if (args.lang != "both" and idioma == "pt") else f"_{idioma}"
        if args.lang == "en":
            sufixo = ""  # Sobrescreve arquivos padrão com textos em inglês se for explicitamente 'en'

        gerar_mcc_tree(destino, lang=idioma, sufixo=sufixo)
        gerar_melhor_zstd(destino, lang=idioma, sufixo=sufixo)
        gerar_figuras_irrc_bzip2(destino, lang=idioma, sufixo=sufixo)

        if args.lang == "both" and idioma == "en":
            # Gera também cópias com sufixo _en para facilitar
            gerar_mcc_tree(destino, lang="en", sufixo="_en")
            gerar_melhor_zstd(destino, lang="en", sufixo="_en")
            gerar_figuras_irrc_bzip2(destino, lang="en", sufixo="_en")

    print("\nOK: Todas as figuras foram geradas com sucesso!")
    print(f"   Faça upload de mcc_tree.png, melhor_zstd.png, irrc_nj.png e bzip2_nj.png no Overleaf.")


if __name__ == "__main__":
    main()

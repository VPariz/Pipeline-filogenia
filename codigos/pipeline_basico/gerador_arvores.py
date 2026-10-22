# ====================================================================
# SCRIPT DE CONSTRUÇÃO DE ÁRVORES FILOGENÉTICAS (NCD-Phylogeny Builder) - V2
# Autor: João Victor Pariz
# Data: 15/11/25 (Modificada para batch processing)
# Propósito: Pipeline automatizado para construir árvores filogenéticas
#            (UPGMA, Linkagem Completa, Neighbor-Joining) a partir de
#            todos os arquivos .txt do diretório, salvando os resultados
#            em uma subpasta dedicada.
# ====================================================================

import os
import sys
import re
import numpy as np
from scipy.cluster.hierarchy import average, complete, to_tree
from scipy.spatial.distance import squareform
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from io import StringIO
from Bio import Phylo

# ====================================================================
# CONSTANTES
# ====================================================================
PASTA_ENTRADA = "resultados_ncd"
PASTA_SAIDA = os.path.join("..", "..", "Resultados", "arvores")

# ====================================================================
# FUNÇÕES DE CONVERSÃO E CLUSTERING (MANTIDAS ORIGINAIS)
# ====================================================================

def converter_para_newick_scipy_ultrametrica(no, nomes_folhas):
    """
    Converte a estrutura da árvore hierárquica (objeto scipy.cluster.hierarchy.to_tree)
    para o formato de string Newick, assumindo natureza ultramétrica (UPGMA/Linkagem Completa).

    Args:
        no: Nó atual da árvore (scipy.cluster.hierarchy.to_tree).
        nomes_folhas (list): Lista ordenada de strings contendo os nomes dos táxons.

    Returns:
        str: Sub-árvore ou folha no formato Newick, incluindo os comprimentos de ramo.
    """
    if no.is_leaf():
        return nomes_folhas[no.id]
    else:
        esquerda, direita = no.get_left(), no.get_right()
        comprimento_esquerda = max((no.dist - esquerda.dist) / 2, 0.0)
        comprimento_direita = max((no.dist - direita.dist) / 2, 0.0)
        return (
            f"({converter_para_newick_scipy_ultrametrica(esquerda, nomes_folhas)}:{comprimento_esquerda:.15f},"
            f"{converter_para_newick_scipy_ultrametrica(direita, nomes_folhas)}:{comprimento_direita:.15f})"
        )

def gerar_newick_upgma(matriz, nomes_folhas):
    """Gera a string Newick para árvore UPGMA (Average Linkage)."""
    matriz_copia = matriz.copy()
    np.fill_diagonal(matriz_copia, 0.0)
    matriz_simetrica = np.minimum(matriz_copia, matriz_copia.T)
    matriz_condensada = squareform(matriz_simetrica, force='tovector')
    Z = average(matriz_condensada)
    raiz = to_tree(Z)
    string_newick = converter_para_newick_scipy_ultrametrica(raiz, nomes_folhas) + ";"
    return string_newick

def gerar_newick_linkagem_completa(matriz, nomes_folhas):
    """Gera a string Newick para árvore de Linkagem Completa (Complete Linkage)."""
    matriz_copia = matriz.copy()
    np.fill_diagonal(matriz_copia, 0.0)
    matriz_simetrica = np.minimum(matriz_copia, matriz_copia.T)
    matriz_condensada = squareform(matriz_simetrica, force='tovector')
    Z = complete(matriz_condensada)
    raiz = to_tree(Z)
    string_newick = converter_para_newick_scipy_ultrametrica(raiz, nomes_folhas) + ";"
    return string_newick

def gerar_newick_neighbor_joining(matriz, nomes_folhas):
    """Gera a string Newick para árvore Neighbor-Joining (NJ) usando Biopython."""
    matriz_copia = matriz.copy()
    n = matriz_copia.shape[0]
    np.fill_diagonal(matriz_copia, 0.0)
    matriz_simetrica = np.minimum(matriz_copia, matriz_copia.T)
    matriz_lower_triangle = [matriz_simetrica[i, :i+1].tolist() for i in range(n)]
    matriz_distancia_bio = DistanceMatrix(nomes_folhas, matriz_lower_triangle)
    construtor = DistanceTreeConstructor()
    arvore_nj = construtor.nj(matriz_distancia_bio)
    handle = StringIO()
    Phylo.write(arvore_nj, handle, "newick")
    string_newick = handle.getvalue().strip()
    return string_newick

# ====================================================================
# FUNÇÃO DE CARREGAMENTO DE MATRIZ NCD (MANTIDA ORIGINAL)
# ====================================================================

def carregar_matriz_ncd(caminho_arquivo):
    """
    Carrega uma matriz de distâncias NCD a partir de um arquivo de texto.
    Retorna (matriz_distancia, nomes_folhas) ou (None, None) em caso de erro.
    """
    padrao = re.compile(r"NCD\((.+?)\.txt,(.+?)\.txt\)=(\S+)")
    distancias = {}
    todos_nomes = set()

    try:
        with open(caminho_arquivo, 'r') as f:
            linhas = f.readlines()
    except FileNotFoundError:
        print(f"Erro: O arquivo '{caminho_arquivo}' não foi encontrado.")
        return None, None

    for linha in linhas:
        match = padrao.search(linha.strip())
        if match:
            nome1 = match.group(1) + '.txt'
            nome2 = match.group(2) + '.txt'
            valor = float(match.group(3))
            distancias[(nome1, nome2)] = valor
            todos_nomes.add(nome1)
            todos_nomes.add(nome2)

    if not todos_nomes:
        print(f"Aviso: Nenhuma distância válida encontrada em '{caminho_arquivo}'. Pulando.")
        return None, None

    nomes_folhas = sorted(list(todos_nomes))
    nome_para_indice = {nome: i for i, nome in enumerate(nomes_folhas)}
    N = len(nomes_folhas)
    matriz = np.zeros((N, N), dtype=float)

    for (nome1, nome2), valor in distancias.items():
        i = nome_para_indice[nome1]
        j = nome_para_indice[nome2]
        matriz[i, j] = valor

    print(f"Carregado {len(distancias)} distâncias para {N} itens.")
    return matriz, nomes_folhas

# ====================================================================
# FUNÇÃO DE SALVAMENTO AUTOMÁTICO (SEM INTERAÇÃO)
# ====================================================================

def salvar_newick_automatico(string_newick, nome_base_arquivo, metodo, pasta_destino):
    """
    Salva a string Newick em um arquivo .nwk dentro da pasta de destino.
    Args:
        string_newick (str): Conteúdo Newick.
        nome_base_arquivo (str): Nome do arquivo original (ex: 'matriz.txt').
        metodo (str): Nome do método (ex: 'upgma', 'complete_linkage', 'neighbor_joining').
        pasta_destino (str): Caminho da pasta onde salvar.
    """
    os.makedirs(pasta_destino, exist_ok=True)
    nome_arquivo_saida = f"{os.path.splitext(nome_base_arquivo)[0]}_{metodo}.nwk"
    caminho_completo = os.path.join(pasta_destino, nome_arquivo_saida)
    with open(caminho_completo, 'w') as f:
        f.write(string_newick)
    print(f"    Salvo: {caminho_completo}")

# ====================================================================
# PROCESSAMENTO EM LOTE PRINCIPAL
# ====================================================================

def main():
    pasta_busca = PASTA_ENTRADA
    if not os.path.exists(pasta_busca) or not [f for f in os.listdir(pasta_busca) if f.endswith('.txt')]:
        alt_pasta = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Resultados", "matrizes"))
        if os.path.exists(alt_pasta) and [f for f in os.listdir(alt_pasta) if f.endswith('.txt')]:
            pasta_busca = alt_pasta
            print(f"Nota: utilizando matrizes NCD pré-computadas de '{alt_pasta}'.\n")
        else:
            print(f"Pasta '{pasta_busca}' não encontrada. Rode o compressor.py primeiro.")
            sys.exit(1)
        
    arquivos_txt = [f for f in os.listdir(pasta_busca) if f.endswith('.txt')]
    if not arquivos_txt:
        print(f"Nenhum arquivo .txt encontrado em '{pasta_busca}'. Encerrando.")
        sys.exit(0)

    # Cria a pasta de saída (se não existir)
    os.makedirs(PASTA_SAIDA, exist_ok=True)
    print(f"Pasta de saída: '{PASTA_SAIDA}' (criada/verificada)\n")

    # Dicionário que mapeia nome do método -> função geradora
    metodos = {
        'upgma': gerar_newick_upgma,
        'complete_linkage': gerar_newick_linkagem_completa,
        'neighbor_joining': gerar_newick_neighbor_joining
    }

    for arquivo in arquivos_txt:
        print(f"\n>>> Processando: {arquivo}")

        # Carrega a matriz NCD
        caminho_arquivo = os.path.join(pasta_busca, arquivo)
        matriz, nomes = carregar_matriz_ncd(caminho_arquivo)
        if matriz is None or nomes is None:
            print(f"  [ERRO] Não foi possível carregar {arquivo}. Pulando.")
            continue

        # Para cada método, gera a árvore e salva automaticamente
        for nome_metodo, func_geradora in metodos.items():
            print(f"  >> Gerando árvore: {nome_metodo.upper()}")
            try:
                newick_str = func_geradora(matriz, nomes)
                if newick_str:
                    salvar_newick_automatico(newick_str, arquivo, nome_metodo, PASTA_SAIDA)
                else:
                    print(f"    [AVISO] {nome_metodo.upper()} retornou string vazia.")
            except Exception as e:
                print(f"    [ERRO] Falha no método {nome_metodo.upper()}: {e}")

        print(f">>> Finalizado: {arquivo}")

    print("\nProcessamento em lote concluído. Todos os arquivos .nwk estão em:", PASTA_SAIDA)

if __name__ == "__main__":
    main()
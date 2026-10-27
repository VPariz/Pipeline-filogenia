# Pipeline de Filogenia via Compressão (NCD)

[![Paper](https://img.shields.io/badge/Paper-BSB_2026-blue?style=flat-square&logo=academia&logoColor=white)](https://github.com/VPariz/Pipeline-filogenia)
[![Poster](https://img.shields.io/badge/Poster-X--Meeting_2026-green?style=flat-square)](https://github.com/VPariz/Pipeline-filogenia)
[![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](LICENSE)

**Projeto de Iniciação Científica (PIBIC/CNPq) — Instituto Militar de Engenharia (IME) & PUC-Rio**  
*Artigo Submetido ao Simpósio Brasileiro de Bioinformática (BSB 2026 / Sociedade Brasileira de Computação - SBC)*

---

## Sobre o Projeto

Este repositório contém o código fonte, dados e pipelines computacionais desenvolvidos para o artigo:  
**"Normalized Compression Distance for Phylogenetic Tree Reconstruction"** (*Pariz, Garcia & Haeusler*).

O trabalho investiga a aplicação da **Distância de Compressão Normalizada (NCD)** — fundamentada na Complexidade de Kolmogorov — para a reconstrução de árvores filogenéticas a partir de genomas mitocondriais completos (mtDNA) de mamíferos sem a necessidade de alinhamento sequencial múltiplo prévio.

### Compressores Avaliados
1. **Compressores Estatísticos Modernos:** `zstd` (nível 22), `brotli` (qualidade 11) e `ppmd` (ordem 16, 64 MB).
2. **Heurística Gramatical Aumentada (ACFG):** Algoritmo `IRR-C` (desenvolvido pelo grupo com identificação de homologia estrutural).
3. **Benchmark Clássico:** `bzip2` (nível -9).

### Algoritmos de Agrupamento Hierárquico
- **UPGMA** (Unweighted Pair Group Method with Arithmetic Mean)
- **Complete Linkage**
- **Neighbor-Joining (NJ)**

### Validação Topológica
As árvores obtidas são validadas estatisticamente contra a **Árvore de Consenso Bayesiana Padrão-Ouro (MCC / VertLife, Upham et al., 2019)** derivada de 10.000 árvores de 31 genes nucleares através das métricas:
- **RF (Robinson-Foulds)** e **nRF (Robinson-Foulds Normalizado)**
- **CTS (Custom Topological Score)** — penalização ponderada de clados ancestrais
- Testes de permutação de **Monte Carlo (p-valor < 0,005)**

---

## Autores & Orientadores

- **João Victor Pariz** (Pesquisador de Iniciação Científica - IME)
- **Prof. Alex de Vasconcellos Garcia, D.Sc.** (Orientador - IME)
- **Prof. Edward Hermann Haeusler, Ph.D.** (Orientador - PUC-Rio)

---

## Estrutura do Repositório

```
Pipeline-filogenia/
├── dados/genomas/            # Sequências genômicas mitocondriais (NCBI GenBank) em .txt
├── Resultados/
│   ├── matrizes/             # Matrizes de distância NCD canônicas pré-computadas (.txt)
│   ├── arvores/              # Árvores Newick (.nwk) para as 15 combinações e árvore MCC
│   └── resumo_metricas.csv   # Tabela consolidada de métricas topológicas (RF, nRF, CTS)
├── codigos/
│   ├── pipeline_basico/
│   │   ├── compressor.py     # Processamento genômico e cálculo das matrizes NCD
│   │   └── gerador_arvores.py# Agrupamento hierárquico e geração de árvores .nwk
│   └── Graficos/
│       ├── main.py           # Interface gráfica de visualização topológica
│       └── gerar_figuras_artigo.py # Script para reproduzir as figuras exatas do paper
└── requirements.txt
```

---

## Como Testar e Reproduzir

### 1. Requisitos
Crie um ambiente virtual com Python 3.9+ e instale as dependências:
```bash
pip install -r requirements.txt
```

### 2. Rodando o Pipeline Completo
Para gerar as matrizes e as árvores do zero:
```bash
cd codigos/pipeline_basico/
python compressor.py
python gerador_arvores.py
```

### 3. Interface Gráfica & Figuras do Artigo
Para abrir a interface gráfica interativa de análise topológica:
```bash
python codigos/Graficos/main.py
```

Para exportar todas as figuras publicadas no artigo:
```bash
python codigos/Graficos/gerar_figuras_artigo.py
```

---

## Citação

Se utilizar este código ou os resultados gerados em sua pesquisa, por favor cite:

```bibtex
@inproceedings{pariz2026ncd,
  title={Normalized Compression Distance for Phylogenetic Tree Reconstruction},
  author={Pariz, Jo{\~a}o Victor and Garcia, Alex de Vasconcellos and Haeusler, Edward Hermann},
  booktitle={Anais do Simp{\'o}sio Brasileiro de Bioinform{\'a}tica (BSB)},
  year={2026},
  publisher={SBC}
}
```

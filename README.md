# NCD Phylogeny Pipeline

**Projeto de Iniciação Científica (PIBIC) — Submetido ao BSB 2026**

Este repositório contém o código fonte e os dados do pipeline computacional desenvolvido para o artigo *"Normalized Compression Distance for Phylogenetic Tree Reconstruction"*, submetido ao evento BSB 2026 (Sociedade Brasileira de Computação). O projeto aplica a **Distância de Compressão Normalizada (NCD)** com compressores padrão (zstd, bzip2, brotli, ppmd) e a heurística gramatical **IRR-C** para inferir filogenias a partir de genomas mitocondriais de mamíferos.

## Autores
- **João Victor Pariz** (Aluno de Iniciação Científica)
- **Prof. Alex de Vasconcellos Garcia** (Orientador - IME)
- **Prof. Edward Hermann Haeusler** (Orientador - PUC-Rio)

---

## Estrutura do Repositório

* `dados/genomas/`: Sequências genômicas mitocondriais originais (mtDNA) em formato texto.
* `Resultados/arvores/`: Árvores filogenéticas pré-computadas (.nwk) e árvore padrão-ouro MCC de referência.
* `codigos/pipeline_basico/`: Scripts computacionais essenciais do pipeline.
  * `compressor.py`: Processa as sequências de DNA e gera as matrizes de distância NCD.
  * `gerador_arvores.py`: Lê as matrizes geradas e constrói as árvores `.nwk` utilizando algoritmos de agrupamento hierárquico (UPGMA, Linkagem Completa, Neighbor-Joining).
* `codigos/Graficos/`: Interface gráfica e módulos de análise topológica. Lê as árvores resultantes, calcula métricas de erro topológico (RF, nRF e CTS) contra o padrão-ouro e plota as imagens e cladogramas utilizados no artigo.

---

## Como testar e reproduzir

### 1. Requisitos
Crie um ambiente virtual com Python 3.9+ e instale as dependências listadas:
```bash
pip install -r requirements.txt
```

As principais bibliotecas utilizadas são: `scipy`, `numpy`, `pandas`, `matplotlib`, `biopython`, `ete3`, e bibliotecas de compressão (`zstandard`, `brotli`, `pyppmd`).

### 2. Rodando o Pipeline Completo
Se desejar gerar as matrizes e as árvores do zero a partir das sequências de genomas, navegue até a pasta `codigos/pipeline_basico/` e execute:
```bash
python compressor.py
python gerador_arvores.py
```
*(Nota: Certifique-se de ter os utilitários de sistema bzip2 instalados para o funcionamento de todas as compressões).*

### 3. Interface de Visualização e Métricas
Para visualizar as árvores geradas, comparar métricas topológicas (RF vs CTS) e reproduzir os gráficos exatos publicados no artigo, utilize a nossa interface visual:
```bash
python codigos/Graficos/main.py
```
Você também pode gerar todas as figuras do artigo automaticamente executando:
```bash
python codigos/Graficos/gerar_figuras_artigo.py
```

---
**Aviso:** Os scripts contêm um modo de operação de diagnóstico/legado (`replica_bat`) focado apenas em rastreabilidade de bugs antigos dos experimentos; **os resultados canônicos validados no artigo derivam do modo `linear_puro`**.

---

## Citação

Se utilizar este código ou os resultados gerados, por favor cite nosso trabalho:

```bibtex
@inproceedings{pariz2026ncd,
  title={Normalized Compression Distance for Phylogenetic Tree Reconstruction},
  author={Pariz, Jo{\~a}o Victor and Garcia, Alex de Vasconcellos and Haeusler, Edward Hermann},
  booktitle={Anais do Simp{\'o}sio Brasileiro de Bioinform{\'a}tica (BSB)},
  year={2026},
  publisher={SBC}
}
```

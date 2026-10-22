# Ferramenta de Análise Filogenética NCD — Versão 4 (PIBIC)

Esta versão evolui a suite de análise filogenética para publicação científica e integração com o novo ecossistema simplificado de ficheiros do projeto PIBIC.

## 📂 Estrutura de Pastas Esperada

```text
Pibic/
├── Resultados/
│   ├── arvores/               <- 20 arquivos Newick limpos
│   ├── experimentais/         <- Experimentos/réplicas
│   ├── matrizes_originais/    <- Matrizes NCD
│   └── relatorio_metricas.txt <- Relatório padronizado
└── codigos/
    └── Graficos/
        ├── main.py
        ├── parser_metricas.py
        ├── visualizador_arvores.py
        ├── comparador.py
        ├── cgi.py             <- Polyfill gerado para Python 3.13
        └── README.md
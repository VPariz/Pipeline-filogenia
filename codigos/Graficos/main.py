# -*- coding: utf-8 -*-
"""
main.py
-------------
GUI principal do nosso projeto PIBIC.
Junta tudo: visualizacao das arvores, comparacao (CTS vs RF)
e geracao de graficos pro artigo.

Pra rodar a interface visual basta executar via terminal:
    python main.py
"""

from __future__ import annotations

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(_THIS_DIR)
sys.path.insert(0, _THIS_DIR)

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import parser_metricas as pm
import visualizador_arvores as va
import comparador as cmp

# --- Pastas padrao (v5) ---
PASTA_PADRAO_ARVORES = os.path.abspath(
    os.path.join(_THIS_DIR, "..", "..", "Resultados", "arvores")
)
PASTA_EXPERIMENTAIS = os.path.abspath(
    os.path.join(_THIS_DIR, "..", "..", "Resultados", "experimentais")
)
PASTA_PADRAO_RELATORIO = os.path.abspath(
    os.path.join(_THIS_DIR, "..", "..", "Resultados")
)
NOME_RELATORIO_PADRAO = "relatorio_metricas.txt"

_NOMES_MODELO_CANDIDATOS = ["modelo.nwk", "artigo.nwk"]


def _initialdir(pasta: str) -> str:
    return pasta if os.path.isdir(pasta) else _THIS_DIR


def _encontrar_arquivo_modelo(pasta: str) -> str | None:
    if not os.path.isdir(pasta):
        return None
    for nome in _NOMES_MODELO_CANDIDATOS:
        caminho = os.path.join(pasta, nome)
        if os.path.exists(caminho):
            return caminho
    for nome in sorted(os.listdir(pasta)):
        if nome.lower().endswith(".nwk") and ("modelo" in nome.lower() or "artigo" in nome.lower()):
            return os.path.join(pasta, nome)
    return None


# ------------------------------------------------------------------------- #
# Paleta de cores
# ------------------------------------------------------------------------- #
COR_BOA = "#d6f5d6"
COR_MEDIA = "#fff6cc"
COR_RUIM = "#f8d7da"
COR_LEGADO = "#e0d4f7"

COR_TEXTO_BOA = "#27ae60"
COR_TEXTO_MEDIA = "#b8860b"
COR_TEXTO_RUIM = "#c0392b"
COR_TEXTO_LEGADO = "#8e44ad"

# Cores para os algoritmos (scatter)
COR_NJ = "#27ae60"       # verde
COR_UPGMA = "#2980b9"    # azul
COR_CL = "#c0392b"       # vermelho
CORES_ALGORITMO = {
    "Neighbor-Joining": COR_NJ,
    "UPGMA": COR_UPGMA,
    "Complete Linkage": COR_CL,
}

# ------------------------------------------------------------------------- #
# Legenda robusta (Bug 1 resolvido)
# ------------------------------------------------------------------------- #
def criar_legenda(parent, itens, **pack_kwargs):
    """
    Monta legenda com blocos solidos '■' coloridos via foreground.
    Nenhum emoji unicode e usado.
    """
    frame = ttk.Frame(parent)
    for cor, texto in itens:
        tk.Label(frame, text="■", foreground=cor, font=("Segoe UI", 10, "bold")).pack(side="left")
        ttk.Label(frame, text=f" {texto}    ", font=("Segoe UI", 8)).pack(side="left")
    pack_kwargs.setdefault("side", "top")
    pack_kwargs.setdefault("fill", "x")
    frame.pack(**pack_kwargs)
    return frame


# ========================================================================= #
# Aba 1 - Visualizador de Arvores (sem alteracoes)
# ========================================================================= #
class AbaVisualizador(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, padding=10)
        self.tree_atual = None
        self.caminho_atual = None
        self._montar_ui()

    def _montar_ui(self):
        barra1 = ttk.Frame(self)
        barra1.pack(side="top", fill="x")
        barra2 = ttk.Frame(self)
        barra2.pack(side="top", fill="x", pady=(4, 8))

        ttk.Button(barra1, text="Abrir arquivo .nwk...", command=self._abrir_arquivo).pack(side="left", padx=(0, 8))

        self.lbl_arquivo = ttk.Label(barra1, text="Nenhum arquivo carregado", font=("Segoe UI", 9, "italic"))
        self.lbl_arquivo.pack(side="left", padx=(0, 16))

        ttk.Button(barra1, text="Exportar (PNG 300dpi + SVG)...", command=self._exportar).pack(side="right")

        self.var_ocultar_outgroup = tk.BooleanVar(value=False)
        ttk.Checkbutton(barra2, text="Ocultar outgroup artificial (randomFasta-N16k)",
                        variable=self.var_ocultar_outgroup, command=self._redesenhar).pack(side="left", padx=(0, 16))

        self.var_usar_distancias = tk.BooleanVar(value=True)
        ttk.Checkbutton(barra2, text="Usar distâncias de ramo (se disponíveis)",
                        variable=self.var_usar_distancias, command=self._redesenhar).pack(side="left", padx=(0, 16))

        self.var_modo_publicacao = tk.BooleanVar(value=False)
        ttk.Checkbutton(barra2, text="Modo publicação (alta qualidade)",
                        variable=self.var_modo_publicacao, command=self._redesenhar).pack(side="left")

        self.fig = Figure(figsize=(9, 7), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)

        self._desenhar_placeholder()

    def _desenhar_placeholder(self):
        self.ax.clear()
        self.ax.text(0.5, 0.5, "Abra um arquivo .nwk para visualizar a árvore",
                     ha="center", va="center", fontsize=12, color="gray")
        self.ax.axis("off")
        self.canvas.draw()

    def _abrir_arquivo(self):
        caminho = filedialog.askopenfilename(
            title="Selecione um arquivo .nwk",
            initialdir=_initialdir(PASTA_PADRAO_ARVORES),
            filetypes=[("Newick tree", "*.nwk"), ("Todos os arquivos", "*.*")],
        )
        if not caminho:
            return
        self.carregar(caminho)

    def carregar(self, caminho: str):
        try:
            self.tree_atual = va.carregar_arvore(caminho)
        except Exception as exc:
            messagebox.showerror("Erro ao ler árvore", f"Não foi possível ler o arquivo:\n{exc}")
            return

        self.caminho_atual = caminho
        nome = os.path.basename(caminho)
        self.lbl_arquivo.config(text=nome)

        if "sem_distancia" in nome.lower():
            self.var_usar_distancias.set(False)

        self._redesenhar()

    def _redesenhar(self):
        if self.tree_atual is None:
            return

        modo_pub = self.var_modo_publicacao.get()
        self.fig.set_size_inches(8, 10) if modo_pub else self.fig.set_size_inches(9, 7)

        self.ax.clear()
        if modo_pub:
            va.desenhar_cladograma(
                self.tree_atual, self.ax,
                titulo=os.path.basename(self.caminho_atual or ""),
                ocultar_outgroup=self.var_ocultar_outgroup.get(),
                usar_distancias=self.var_usar_distancias.get(),
                modo_publicacao=True,
            )
        else:
            va.desenhar_cladograma(
                self.tree_atual, self.ax,
                titulo=os.path.basename(self.caminho_atual or ""),
                ocultar_outgroup=self.var_ocultar_outgroup.get(),
                usar_distancias=self.var_usar_distancias.get(),
            )
        self.fig.tight_layout()
        self.canvas.draw()

    def _exportar(self):
        if self.tree_atual is None:
            messagebox.showinfo("Nada para exportar", "Carregue uma árvore primeiro.")
            return
        destino = filedialog.asksaveasfilename(
            title="Exportar árvore (será salvo como .png e .svg)",
            defaultextension=".png",
            filetypes=[("Imagem PNG", "*.png")],
            initialfile=os.path.splitext(os.path.basename(self.caminho_atual))[0],
            initialdir=_initialdir(PASTA_PADRAO_ARVORES),
        )
        if not destino:
            return
        base_sem_ext = os.path.splitext(destino)[0]
        va.exportar_figura(self.fig, base_sem_ext)
        messagebox.showinfo("Exportado", f"Arquivos salvos:\n{base_sem_ext}.png\n{base_sem_ext}.svg")


# ========================================================================= #
# Aba 2 - Analisador de Relatorio de Metricas (v5 com CTS)
# ========================================================================= #
class AbaMetricas(ttk.Frame):
    COLUNAS = [
        ("arquivo", "Arquivo", 240),
        ("compressor", "Compressor", 100),
        ("abordagem", "Abordagem", 110),
        ("algoritmo", "Algoritmo", 130),
        ("nrf", "nRF", 80),
        ("cts", "CTS", 90),           # Antiga "diff_prof"
        ("rf_absoluto", "RF Absoluto", 90),
        ("pearson_r", "Pearson r", 90),
        ("spearman_r", "Spearman r", 90),
        ("erro_escala", "Erro Escala", 90),
        ("p_value_topologia", "p-value (Topologia)", 130),
        ("p_value_distancias", "p-value (Distâncias)", 130),
    ]

    def __init__(self, master, callback_abrir_no_visualizador=None):
        super().__init__(master, padding=10)
        self.df_completo = None
        self.callback_abrir_no_visualizador = callback_abrir_no_visualizador
        self._sort_col = None
        self._sort_desc = False
        self._caminho_relatorio_estatico = None
        self.modo_dinamico = False
        self.caminho_referencia_atual = None
        self._montar_ui()
        self._tentar_autocarregar()

    def _montar_ui(self):
        # Barra de referencia dinamica (v3)
        barra_ref = ttk.Frame(self)
        barra_ref.pack(side="top", fill="x", pady=(0, 6))

        ttk.Label(barra_ref, text="Árvore de referência:", font=("Segoe UI", 9, "bold")).pack(side="left")
        self.lbl_referencia_atual = ttk.Label(barra_ref, text="Modelo Taxonômico (relatório estático)")
        self.lbl_referencia_atual.pack(side="left", padx=(4, 14))

        ttk.Button(barra_ref, text="Mudar Referência...", command=self._mudar_referencia).pack(side="left", padx=(0, 6))
        ttk.Button(barra_ref, text="Restaurar Padrão", command=self._restaurar_padrao).pack(side="left", padx=(0, 14))

        self.var_calcular_pvalue = tk.BooleanVar(value=False)
        ttk.Checkbutton(barra_ref, text="Calcular significância estatística (p-value)",
                        variable=self.var_calcular_pvalue).pack(side="left")

        self.progress = ttk.Progressbar(barra_ref, mode="determinate", length=200)

        # Filtros
        barra = ttk.Frame(self)
        barra.pack(side="top", fill="x", pady=(0, 8))

        ttk.Button(barra, text="Abrir relatório de métricas...", command=self._abrir_relatorio).pack(side="left", padx=(0, 12))
        ttk.Button(barra, text="Salvar Relatório Dinâmico...", command=self._salvar_relatorio).pack(side="left", padx=(0, 12))

        ttk.Label(barra, text="Compressor:").pack(side="left")
        self.cb_compressor = ttk.Combobox(barra, state="readonly", width=14, values=["Todos"])
        self.cb_compressor.set("Todos")
        self.cb_compressor.pack(side="left", padx=(2, 10))
        self.cb_compressor.bind("<<ComboboxSelected>>", lambda e: self._aplicar_filtros())

        ttk.Label(barra, text="Algoritmo:").pack(side="left")
        self.cb_algoritmo = ttk.Combobox(barra, state="readonly", width=16, values=["Todos"])
        self.cb_algoritmo.set("Todos")
        self.cb_algoritmo.pack(side="left", padx=(2, 10))
        self.cb_algoritmo.bind("<<ComboboxSelected>>", lambda e: self._aplicar_filtros())

        ttk.Label(barra, text="Abordagem:").pack(side="left")
        self.cb_abordagem = ttk.Combobox(barra, state="readonly", width=14, values=["Todas"])
        self.cb_abordagem.set("Todas")
        self.cb_abordagem.pack(side="left", padx=(2, 10))
        self.cb_abordagem.bind("<<ComboboxSelected>>", lambda e: self._aplicar_filtros())

        ttk.Label(barra, text="Top:").pack(side="left")
        self.cb_top = ttk.Combobox(barra, state="readonly", width=6, values=["Todos", "3", "5", "10"])
        self.cb_top.set("Todos")
        self.cb_top.pack(side="left", padx=(2, 8))
        self.cb_top.bind("<<ComboboxSelected>>", lambda e: self._aplicar_filtros())

        ttk.Button(barra, text="Restaurar ordenação padrão (CTS + nRF)",
                   command=self._restaurar_ordenacao_padrao).pack(side="left")

        # Tabela
        cols = [c[0] for c in self.COLUNAS]
        self.tree = ttk.Treeview(self, columns=cols, show="headings", selectmode="browse")
        for col_id, titulo, largura in self.COLUNAS:
            self.tree.heading(col_id, text=titulo, command=lambda c=col_id: self._ordenar_por(c))
            self.tree.column(col_id, width=largura, anchor="center" if col_id != "arquivo" else "w")

        vsb = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side="top", fill="both", expand=True)
        vsb.pack(side="right", fill="y")

        self.tree.tag_configure("boa", background=COR_BOA)
        self.tree.tag_configure("media", background=COR_MEDIA)
        self.tree.tag_configure("ruim", background=COR_RUIM)
        self.tree.tag_configure("legado", background=COR_LEGADO)

        self.tree.bind("<Double-1>", self._duplo_clique)

        # Rodape com legenda (v5: CTS)
        rodape = ttk.Frame(self)
        rodape.pack(side="bottom", fill="x", pady=(6, 0))
        self.lbl_status = ttk.Label(rodape, text="Nenhum relatório carregado.")
        self.lbl_status.pack(side="left")

        criar_legenda(
            rodape,
            [
                (COR_TEXTO_BOA, "CTS ≤ 3.0 e nRF ≤ 0.12"),
                (COR_TEXTO_MEDIA, "CTS ≤ 5.0 ou nRF ≤ 0.18"),
                (COR_TEXTO_RUIM, "demais"),
                (COR_TEXTO_LEGADO, "experimentais (LEGADO)"),
            ],
            side="right",
        )

        ttk.Label(self, text="Duplo-clique numa linha para ver a árvore correspondente.",
                  font=("Segoe UI", 8), foreground="#666666").pack(side="bottom", anchor="w", pady=(2, 0))

        self.lbl_pastas = ttk.Label(
            self,
            text=f"Pasta padrão de árvores: {PASTA_PADRAO_ARVORES}\n"
                 f"Pasta de experimentais: {PASTA_EXPERIMENTAIS}\n"
                 f"Pasta padrão de relatórios: {PASTA_PADRAO_RELATORIO}",
            font=("Segoe UI", 8), foreground="#555555", justify="left",
        )
        self.lbl_pastas.pack(side="bottom", fill="x", anchor="w", pady=(4, 0))

    def _tentar_autocarregar(self):
        caminho = os.path.join(PASTA_PADRAO_RELATORIO, NOME_RELATORIO_PADRAO)
        if os.path.exists(caminho):
            try:
                self._carregar(caminho)
                self._caminho_relatorio_estatico = caminho
                self.lbl_status.config(text=f"Relatório carregado automaticamente: {caminho}")
            except Exception:
                pass

    def _abrir_relatorio(self):
        caminho = filedialog.askopenfilename(
            title="Selecione o relatório de métricas (.txt)",
            initialdir=_initialdir(PASTA_PADRAO_RELATORIO),
            filetypes=[("Relatório de texto", "*.txt"), ("Todos os arquivos", "*.*")],
        )
        if not caminho:
            return
        self._carregar(caminho)
        self._caminho_relatorio_estatico = caminho
        self.modo_dinamico = False
        self.caminho_referencia_atual = None
        self.lbl_referencia_atual.config(text="Modelo Taxonômico (relatório estático)")

    def _carregar(self, caminho):
        df = pm.parse_relatorio(caminho)
        if df.empty:
            messagebox.showwarning("Relatório vazio", "Nenhum bloco 'ARQUIVO:' foi encontrado.")
            return
        self.df_completo = df
        self._popular_combos()
        self._sort_col = None
        self._aplicar_filtros()

    # --- v3/v5: referencia dinamica ---
    def _mudar_referencia(self):
        caminho = filedialog.askopenfilename(
            title="Selecione a árvore de referência para recalcular tudo dinamicamente",
            initialdir=_initialdir(PASTA_PADRAO_ARVORES),
            filetypes=[("Newick tree", "*.nwk"), ("Todos os arquivos", "*.*")],
        )
        if not caminho:
            return
        self._recalcular_dinamico(caminho)

    def _restaurar_padrao(self):
        if self._caminho_relatorio_estatico and os.path.exists(self._caminho_relatorio_estatico):
            self._carregar(self._caminho_relatorio_estatico)
            self.modo_dinamico = False
            self.caminho_referencia_atual = None
            self.lbl_referencia_atual.config(text="Modelo Taxonômico (relatório estático)")
            self.lbl_status.config(text=f"Relatório estático restaurado: {self._caminho_relatorio_estatico}")
        else:
            messagebox.showinfo("Relatório padrão indisponível", "Nenhum relatório estático foi carregado.")

    def _recalcular_dinamico(self, caminho_ref):
        if not os.path.isdir(PASTA_PADRAO_ARVORES):
            messagebox.showerror("Pasta de árvores não encontrada", f"A pasta {PASTA_PADRAO_ARVORES} não existe.")
            return

        self.progress.pack(side="left", padx=(10, 0))
        self.progress["value"] = 0

        def callback(i, total, nome):
            self.progress["maximum"] = max(total, 1)
            self.progress["value"] = i
            self.lbl_status.config(text=f"Calculando dinamicamente {i}/{total}: {nome}")
            self.update_idletasks()

        try:
            df = pm.calcular_relatorio_dinamico(
                caminho_ref, PASTA_PADRAO_ARVORES,
                calcular_pvalue=self.var_calcular_pvalue.get(),
                n_baseline=30,
                progress_callback=callback,
            )
        except Exception as exc:
            messagebox.showerror("Erro no cálculo dinâmico", str(exc))
            self.progress.pack_forget()
            return

        self.progress.pack_forget()

        if df.empty:
            messagebox.showwarning("Nenhum resultado", "Nenhuma árvore pôde ser comparada contra essa referência.")
            return

        n_erros = int(df["erro_calculo"].notna().sum()) if "erro_calculo" in df.columns else 0

        self.df_completo = df
        self.modo_dinamico = True
        self.caminho_referencia_atual = caminho_ref
        self.lbl_referencia_atual.config(
            text=f"{os.path.basename(caminho_ref)}  (dinâmico — aproximado, ver README)"
        )
        self._popular_combos()
        self._sort_col = None
        self._aplicar_filtros()

        msg = f"Cálculo dinâmico concluído: {len(df)} árvore(s) comparada(s)."
        if n_erros:
            msg += f" ({n_erros} com erro de cálculo)"
        self.lbl_status.config(text=msg)

    # --- filtros e ordenacao ---
    def _popular_combos(self):
        df = self.df_completo
        self.cb_compressor["values"] = ["Todos"] + sorted(df["compressor"].dropna().unique().tolist())
        self.cb_algoritmo["values"] = ["Todos"] + sorted(df["algoritmo"].dropna().unique().tolist())
        abordagens = sorted(a for a in df["abordagem"].dropna().unique().tolist() if a)
        self.cb_abordagem["values"] = ["Todas"] + abordagens

    def _aplicar_filtros(self):
        if self.df_completo is None:
            return
        df = self.df_completo.copy()

        if self.cb_compressor.get() != "Todos":
            df = df[df["compressor"] == self.cb_compressor.get()]
        if self.cb_algoritmo.get() != "Todos":
            df = df[df["algoritmo"] == self.cb_algoritmo.get()]
        if self.cb_abordagem.get() != "Todas":
            df = df[df["abordagem"] == self.cb_abordagem.get()]

        df = self._ordenar_df(df)

        top_sel = self.cb_top.get()
        if top_sel != "Todos":
            df = df.sort_values(["cts", "nrf"], ascending=[True, True], na_position="last")
            df = df.head(int(top_sel))
            df = self._ordenar_df(df)

        self._popular_tabela(df)

    def _ordenar_df(self, df):
        if self._sort_col is None:
            return df.sort_values(["cts", "nrf"], ascending=[True, True], na_position="last")
        return df.sort_values(self._sort_col, ascending=not self._sort_desc, na_position="last")

    def _ordenar_por(self, col_id):
        if self._sort_col == col_id:
            self._sort_desc = not self._sort_desc
        else:
            self._sort_col = col_id
            self._sort_desc = False
        self._aplicar_filtros()

    def _restaurar_ordenacao_padrao(self):
        self._sort_col = None
        self._sort_desc = False
        self._aplicar_filtros()

    def _popular_tabela(self, df):
        self.tree.delete(*self.tree.get_children())
        for _, row in df.iterrows():
            valores = []
            for col_id, _, _ in self.COLUNAS:
                v = row.get(col_id)
                if isinstance(v, float):
                    v = f"{v:.4f}"
                valores.append(v if v is not None else "—")

            tag = pm.classificar_linha(row)   # usa "cts"
            self.tree.insert("", "end", iid=row["arquivo"], values=valores, tags=(tag,))

        self.lbl_status.config(text=f"{len(df)} arquivo(s) exibido(s) de {len(self.df_completo)} no relatório.")

    def _duplo_clique(self, event):
        item = self.tree.identify_row(event.y)
        if not item or self.callback_abrir_no_visualizador is None:
            return
        self.callback_abrir_no_visualizador(item)

    def _salvar_relatorio(self):
        if self.df_completo is None:
            messagebox.showinfo("Nada para salvar", "Carregue ou calcule um relatório primeiro.")
            return

        destino = filedialog.asksaveasfilename(
            title="Salvar relatório (texto ou Excel)",
            defaultextension=".txt",
            filetypes=[("Texto", "*.txt"), ("Planilha Excel", "*.xlsx")],
            initialdir=_initialdir(PASTA_PADRAO_RELATORIO),
            initialfile="relatorio_dinamico" if self.modo_dinamico else "relatorio_exportado",
        )
        if not destino:
            return

        referencia = (
            os.path.basename(self.caminho_referencia_atual)
            if self.caminho_referencia_atual else "Modelo Taxonômico (relatório estático)"
        )

        try:
            if destino.lower().endswith(".xlsx"):
                self.df_completo.to_excel(destino, index=False)
            else:
                texto = pm.formatar_relatorio_texto(self.df_completo, referencia=referencia)
                with open(destino, "w", encoding="utf-8") as f:
                    f.write(texto)
        except ImportError:
            messagebox.showerror("Dependência ausente", "Para exportar .xlsx instale openpyxl")
            return
        except Exception as exc:
            messagebox.showerror("Erro ao salvar", str(exc))
            return

        messagebox.showinfo("Salvo", f"Relatório salvo em:\n{destino}")


# ========================================================================= #
# Aba 3 - Comparador de Arvores (Bug 2 corrigido; v5 com CTS)
# ========================================================================= #
class AbaComparador(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, padding=10)
        self.caminho_a = None
        self.caminho_b = None
        self.tree_a = None
        self.tree_b = None
        self.analise_atual = None
        self.calculo_completo_atual = None
        self._montar_ui()
        self._tentar_pre_selecionar_modelo()

    def _montar_ui(self):
        barra = ttk.Frame(self)
        barra.pack(side="top", fill="x", pady=(0, 6))

        col_a = ttk.Frame(barra)
        col_a.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ttk.Label(col_a, text="Árvore A (referência — sempre preta)", font=("Segoe UI", 9, "bold")).pack(anchor="w")
        linha_a = ttk.Frame(col_a)
        linha_a.pack(fill="x")
        ttk.Button(linha_a, text="Selecionar...", command=lambda: self._selecionar("a")).pack(side="left")
        self.lbl_a = ttk.Label(linha_a, text="(nenhum arquivo)")
        self.lbl_a.pack(side="left", padx=8)

        col_b = ttk.Frame(barra)
        col_b.pack(side="left", fill="x", expand=True, padx=(10, 0))
        ttk.Label(col_b, text="Árvore B (resultado a comparar — colorida)", font=("Segoe UI", 9, "bold")).pack(anchor="w")
        linha_b = ttk.Frame(col_b)
        linha_b.pack(fill="x")
        ttk.Button(linha_b, text="Selecionar...", command=lambda: self._selecionar("b")).pack(side="left")
        self.lbl_b = ttk.Label(linha_b, text="(nenhum arquivo)")
        self.lbl_b.pack(side="left", padx=8)

        barra_botoes = ttk.Frame(self)
        barra_botoes.pack(side="top", fill="x", pady=(4, 6))
        ttk.Button(barra_botoes, text="Comparar (RF - bipartições)",
                   command=lambda: self._comparar(modo="rf")).pack(side="left", padx=(0, 8))
        ttk.Button(barra_botoes, text="Comparar (CTS - clados)",
                   command=lambda: self._comparar(modo="cts")).pack(side="left", padx=(0, 8))

        self.var_calcular_pvalue = tk.BooleanVar(value=False)
        ttk.Checkbutton(barra_botoes, text="Calcular significância estatística (p-value)",
                        variable=self.var_calcular_pvalue).pack(side="left", padx=(8, 0))

        ttk.Button(barra_botoes, text="Exportar comparação (SVG+PNG)...",
                   command=self._exportar_comparacao).pack(side="right")

        self.lbl_metricas = ttk.Label(self, text="", font=("Segoe UI", 10), justify="left")
        self.lbl_metricas.pack(side="top", fill="x", pady=(0, 4))

        area_figuras = ttk.Frame(self)
        area_figuras.pack(side="top", fill="both", expand=True)

        self.fig_a = Figure(figsize=(5.5, 7), dpi=100)
        self.ax_a = self.fig_a.add_subplot(111)
        self.canvas_a = FigureCanvasTkAgg(self.fig_a, master=area_figuras)
        self.canvas_a.get_tk_widget().pack(side="left", fill="both", expand=True)

        self.fig_b = Figure(figsize=(5.5, 7), dpi=100)
        self.ax_b = self.fig_b.add_subplot(111)
        self.canvas_b = FigureCanvasTkAgg(self.fig_b, master=area_figuras)
        self.canvas_b.get_tk_widget().pack(side="left", fill="both", expand=True)

        # Legenda Tkinter (v5: sem Diff Prof)
        criar_legenda(
            self,
            [
                ("#2ecc71", "Clado correto (presente no modelo)"),
                ("#e74c3c", "Clado errado (ausente no modelo)"),
            ],
            side="top", pady=(6, 0),
        )

        self._placeholder()

    def _placeholder(self):
        for ax, canvas, letra in ((self.ax_a, self.canvas_a, "A"), (self.ax_b, self.canvas_b, "B")):
            ax.clear()
            ax.text(0.5, 0.5, f"Selecione a árvore {letra}", ha="center", va="center", color="gray")
            ax.axis("off")
            canvas.draw()

    def _tentar_pre_selecionar_modelo(self):
        candidato = _encontrar_arquivo_modelo(PASTA_PADRAO_ARVORES)
        if candidato:
            self._carregar_slot("a", candidato)

    def _selecionar(self, slot):
        caminho = filedialog.askopenfilename(
            title=f"Selecione a árvore {slot.upper()}",
            initialdir=_initialdir(PASTA_PADRAO_ARVORES),
            filetypes=[("Newick tree", "*.nwk"), ("Todos os arquivos", "*.*")],
        )
        if not caminho:
            return
        self._carregar_slot(slot, caminho)

    def _carregar_slot(self, slot, caminho):
        if slot == "a":
            self.caminho_a = caminho
            self.lbl_a.config(text=os.path.basename(caminho))
        else:
            self.caminho_b = caminho
            self.lbl_b.config(text=os.path.basename(caminho))

        # Atualiza a comparacao e redesenha a tela imediatamente se ambas estiverem selecionadas
        if self.caminho_a and self.caminho_b:
            self._comparar(modo="cts")


    def _comparar(self, modo: str):
        if not self.caminho_a or not self.caminho_b:
            messagebox.showinfo("Selecione as árvores", "Escolha os arquivos das árvores A e B primeiro.")
            return

        try:
            self.tree_a = va.carregar_arvore(self.caminho_a)
            self.tree_b = va.carregar_arvore(self.caminho_b)
        except Exception as exc:
            messagebox.showerror("Erro ao ler árvores", str(exc))
            return

        try:
            self.analise_atual = cmp.analisar_biparticoes(self.caminho_a, self.caminho_b)
        except Exception as exc:
            messagebox.showerror("Erro ao comparar topologias", str(exc))
            return

        self.calculo_completo_atual = cmp.executar_comparacao_completa(
            self.caminho_a, self.caminho_b,
            calcular_pvalue=self.var_calcular_pvalue.get(), n_baseline=30,
        )

        # Arvore A: preta
        va.desenhar_cladograma(self.tree_a, self.ax_a, titulo=os.path.basename(self.caminho_a),
                               usar_distancias=False)
        self.fig_a.tight_layout()
        self.canvas_a.draw()

        # Arvore B: colorida; contagem no titulo (Bug 2)
        if modo == "rf":
            color_fn = cmp.gerador_cor_rf(self.analise_atual)
            corretos, total = self.analise_atual["corretos_rf"], self.analise_atual["total_b"]
            titulo_extra = f"{corretos}/{total} clados corretos (RF)"
        else:  # modo == "cts"
            color_fn = cmp.gerador_cor_cts(self.analise_atual)   # antigo diffprof
            corretos, total = self.analise_atual["corretos_cts"], self.analise_atual["total_b"]
            titulo_extra = f"{corretos}/{total} clados corretos (CTS)"

        va.desenhar_cladograma(self.tree_b, self.ax_b, titulo=os.path.basename(self.caminho_b),
                               usar_distancias=False, color_fn=color_fn,
                               titulo_extra=titulo_extra)
        self.fig_b.tight_layout()
        self.canvas_b.draw()

        self._atualizar_painel_metricas(modo)

    def _atualizar_painel_metricas(self, modo):
        a = self.analise_atual
        c = self.calculo_completo_atual

        if a["erro_rf"]:
            texto = f"⚠ Não foi possível calcular RF/nRF: {a['erro_rf']}"
        else:
            texto = (
                f"RF Absoluto: {a['rf_absoluto']:.0f}  |  nRF: {a['nrf']:.4f}  |  "
                f"Clados comuns (RF): {a['corretos_rf']}/{a['total_b']}  |  "
                f"Clados corretos (CTS): {a['corretos_cts']}/{a['total_b']}"
            )
            if a["taxons_exclusivos_a"] or a["taxons_exclusivos_b"]:
                texto += f"\nExclusivos de A: {a['taxons_exclusivos_a']}   Exclusivos de B: {a['taxons_exclusivos_b']}"
            texto += f"\n(Visão exibida na árvore B: {'RF/bipartições' if modo == 'rf' else 'CTS/clados enraizados'})"

        if c and not c.get("erro"):
            texto += (
                f"\n\nCálculo dinâmico numérico (aproximado — ver README):  "
                f"CTS: {c['cts']:.4f}  |  Pearson r: {c['pearson_r']:.4f}  |  "
                f"Spearman r: {c['spearman_r']:.4f}  |  Erro Escala: {c['erro_escala']:.4f}"
            )
            if c.get("p_value_topologia") is not None:
                texto += f"\np-value (Topologia): {c['p_value_topologia']:.4f}   " \
                         f"p-value (Distâncias): {c['p_value_distancias']:.4f}"
        elif c and c.get("erro"):
            texto += f"\n\n⚠ Cálculo numérico dinâmico falhou: {c['erro']}"

        self.lbl_metricas.config(text=texto)

    def _exportar_comparacao(self):
        if self.analise_atual is None:
            messagebox.showinfo("Nada para exportar", "Faça uma comparação primeiro.")
            return
        destino = filedialog.asksaveasfilename(
            title="Exportar comparação (serão salvos _A e _B em .png e .svg)",
            defaultextension=".png",
            filetypes=[("Imagem PNG", "*.png")],
            initialfile="comparacao",
            initialdir=_initialdir(PASTA_PADRAO_ARVORES),
        )
        if not destino:
            return
        base = os.path.splitext(destino)[0]
        va.exportar_figura(self.fig_a, base + "_A")
        va.exportar_figura(self.fig_b, base + "_B")
        messagebox.showinfo("Exportado", f"Arquivos salvos:\n{base}_A.png / .svg\n{base}_B.png / .svg")


# ========================================================================= #
# Aba 4 - Graficos Comparativos (v5 com CTS)
# ========================================================================= #
class AbaGraficos(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, padding=10)
        self.df_atual = None
        self.metrica_atual = "cts"   # padrao
        self.incluir_experimentais = False
        self._montar_ui()
        self._carregar_dados_padrao()

    def _montar_ui(self):
        barra = ttk.Frame(self)
        barra.pack(side="top", fill="x", pady=(0, 8))

        ttk.Label(barra, text="Métrica:").pack(side="left")
        self.cb_metrica = ttk.Combobox(
            barra, state="readonly", width=24,
            values=["CTS (Custom Topological Score)", "nRF", "Pearson r", "Spearman r", "Erro Escala"]
        )
        self.cb_metrica.set("CTS (Custom Topological Score)")
        self.cb_metrica.pack(side="left", padx=(2, 12))
        self.cb_metrica.bind("<<ComboboxSelected>>", lambda e: self._atualizar_graficos())

        self.var_incluir_exp = tk.BooleanVar(value=False)
        ttk.Checkbutton(barra, text="Incluir experimentos adicionais",
                        variable=self.var_incluir_exp, command=self._atualizar_graficos).pack(side="left", padx=(0, 12))

        self.var_lang_en = tk.BooleanVar(value=False)
        ttk.Checkbutton(barra, text="Gráficos em Inglês (English)",
                        variable=self.var_lang_en, command=self._atualizar_graficos).pack(side="left", padx=(0, 12))

        ttk.Button(barra, text="Recarregar dados", command=self._carregar_dados_padrao).pack(side="left", padx=(0, 12))
        ttk.Button(barra, text="Exportar Gráficos (PNG 300dpi + SVG)...",
                   command=self._exportar_graficos).pack(side="left", padx=(0, 12))
        ttk.Button(barra, text="Exportar Tabela (Excel)...",
                   command=self._exportar_tabela).pack(side="left")

        self.lbl_status = ttk.Label(self, text="", font=("Segoe UI", 8), foreground="#555555")
        self.lbl_status.pack(side="top", fill="x", pady=(0, 4))

        # Area dos graficos: dois paineis lado a lado
        area_plots = ttk.Frame(self)
        area_plots.pack(side="top", fill="both", expand=True)

        # Scatter plot (esquerda)
        self.fig_scatter = Figure(figsize=(6, 6), dpi=100)
        self.ax_scatter = self.fig_scatter.add_subplot(111)
        self.canvas_scatter = FigureCanvasTkAgg(self.fig_scatter, master=area_plots)
        self.canvas_scatter.get_tk_widget().pack(side="left", fill="both", expand=True)

        # Bar chart (direita)
        self.fig_bar = Figure(figsize=(6, 6), dpi=100)
        self.ax_bar = self.fig_bar.add_subplot(111)
        self.canvas_bar = FigureCanvasTkAgg(self.fig_bar, master=area_plots)
        self.canvas_bar.get_tk_widget().pack(side="left", fill="both", expand=True)

        # Legenda de cores (Bug 1: sem emojis)
        criar_legenda(
            self,
            [
                (COR_NJ, "Neighbor-Joining"),
                (COR_UPGMA, "UPGMA"),
                (COR_CL, "Complete Linkage"),
                ("#555555", "○ compressores padrão (bzip2, zstd, brotli, ppmd)"),
                ("#555555", "★ Ficd/21g e IRR-C/aFicd"),
            ],
            side="bottom", pady=(4, 0),
        )

        # Texto de metadados
        ttk.Label(self, text="Scatter: nRF (X) vs métrica selecionada (Y) · Barras: métrica agrupada por compressor",
                  font=("Segoe UI", 8), foreground="#666666").pack(side="bottom", anchor="w", pady=(2, 0))

    def _carregar_dados_padrao(self):
        caminho_rel = os.path.join(PASTA_PADRAO_RELATORIO, NOME_RELATORIO_PADRAO)
        if os.path.exists(caminho_rel):
            try:
                self.df_atual = pm.parse_relatorio(caminho_rel)
                self.lbl_status.config(text=f"Dados carregados de {caminho_rel}")
                self._atualizar_graficos()
                return
            except Exception as e:
                self.lbl_status.config(text=f"Erro ao carregar relatório: {e}")

        modelo = _encontrar_arquivo_modelo(PASTA_PADRAO_ARVORES)
        if modelo:
            try:
                self.lbl_status.config(text="Recalculando dinamicamente...")
                self.update_idletasks()
                self.df_atual = pm.calcular_relatorio_dinamico(
                    modelo, PASTA_PADRAO_ARVORES,
                    calcular_pvalue=False, n_baseline=30,
                    progress_callback=None,
                )
                self.lbl_status.config(text=f"Dados dinâmicos calculados (referência: {os.path.basename(modelo)})")
                self._atualizar_graficos()
                return
            except Exception as e:
                self.lbl_status.config(text=f"Erro no cálculo dinâmico: {e}")

        self.df_atual = None
        self.lbl_status.config(text="⚠ Nenhum dado disponível. Verifique a pasta de árvores.")
        self._placeholder_plots()

    def _placeholder_plots(self):
        for ax, fig in ((self.ax_scatter, self.fig_scatter), (self.ax_bar, self.fig_bar)):
            ax.clear()
            ax.text(0.5, 0.5, "Sem dados para exibir", ha="center", va="center", color="gray", fontsize=14)
            ax.axis("off")
            fig.tight_layout()
        self.canvas_scatter.draw()
        self.canvas_bar.draw()

    def _atualizar_graficos(self):
        if self.df_atual is None or self.df_atual.empty:
            self._placeholder_plots()
            return

        # Mapeia nome da metrica para coluna (v5: CTS)
        mapa_metrica = {
            "CTS (Custom Topological Score)": "cts",
            "nRF": "nrf",
            "Pearson r": "pearson_r",
            "Spearman r": "spearman_r",
            "Erro Escala": "erro_escala",
        }
        col_metrica = mapa_metrica.get(self.cb_metrica.get(), "cts")
        self.metrica_atual = col_metrica

        df = self.df_atual.copy()

        # Filtra linhas com dados validos para a metrica e nrf
        df = df[df[col_metrica].notna() & df["nrf"].notna()]

        # Separa principais vs experimentais
        df_princ = df[~df.get("legado", False)]
        if self.var_incluir_exp.get():
            df_exp = df[df.get("legado", False)]
        else:
            df_exp = pd.DataFrame()

        # --- Scatter plot: nRF (X) vs metrica (Y) ---
        self.ax_scatter.clear()
        self._desenhar_scatter(self.ax_scatter, df_princ, df_exp, col_metrica)
        self.fig_scatter.tight_layout()
        self.canvas_scatter.draw()

        # --- Bar chart: metrica por compressor, agrupado por algoritmo ---
        self.ax_bar.clear()
        self._desenhar_barras(self.ax_bar, df_princ, col_metrica)
        self.fig_bar.tight_layout()
        self.canvas_bar.draw()

    def _desenhar_scatter(self, ax, df_princ, df_exp, col_metrica):
        is_en = getattr(self, "var_lang_en", None) and self.var_lang_en.get()
        lbl_x = "nRF (lower = better)" if is_en else "nRF (menor = melhor)"
        lbl_better = " (lower = better)" if is_en else " (menor = melhor)"
        lbl_title = "Metrics Comparison: nRF vs " if is_en else "Comparação de Métricas: nRF vs "

        ax.set_xlabel(lbl_x, fontsize=10)
        if col_metrica == "cts":
            ax.set_ylabel("CTS (Custom Topological Score)" + lbl_better, fontsize=10)
        else:
            ax.set_ylabel(self.cb_metrica.get() + lbl_better, fontsize=10)
        ax.set_title(lbl_title + self.cb_metrica.get(), fontsize=11, fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.6)

        for _, row in df_princ.iterrows():
            compressor = row.get("compressor", "")
            algoritmo = row.get("algoritmo", "")
            x = row["nrf"]
            y = row[col_metrica]

            cor = CORES_ALGORITMO.get(algoritmo, "#333333")
            is_ficd = compressor in ("Ficd/21g", "IRR-C/aFicd")
            marker = "*" if is_ficd else "o"
            size = 120 if is_ficd else 80

            ax.scatter(x, y, c=cor, marker=marker, s=size, edgecolors="black", linewidth=0.5, alpha=0.9)

            rotulo = self._abreviatura_compressor(compressor)
            ax.annotate(rotulo, (x, y), xytext=(4, 4), textcoords="offset points",
                        fontsize=7, alpha=0.8, fontweight="bold")

        if not df_exp.empty:
            lbl_exp_legend = "Experimental (LEGACY)" if (getattr(self, "var_lang_en", None) and self.var_lang_en.get()) else "Experimentais (LEGADO)"
            ax.scatter(df_exp["nrf"], df_exp[col_metrica],
                       c="#8e44ad", marker="s", s=60, edgecolors="black", linewidth=0.5,
                       alpha=0.6, label=lbl_exp_legend)
            for _, row in df_exp.iterrows():
                rotulo = self._abreviatura_compressor(row.get("compressor", ""))
                ax.annotate(rotulo, (row["nrf"], row[col_metrica]), xytext=(4, -8),
                            textcoords="offset points", fontsize=6, alpha=0.6, color="#8e44ad")

        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=COR_NJ, label="Neighbor-Joining"),
            Patch(facecolor=COR_UPGMA, label="UPGMA"),
            Patch(facecolor=COR_CL, label="Complete Linkage"),
        ]
        ax.legend(handles=legend_elements, loc="best", fontsize=8, framealpha=0.8)

        x_vals = df_princ["nrf"].dropna()
        y_vals = df_princ[col_metrica].dropna()
        if not x_vals.empty and not y_vals.empty:
            x_max = x_vals.max() * 1.15
            y_max = y_vals.max() * 1.15
            ax.set_xlim(-0.02, max(0.5, x_max) if x_max > 0 else 0.5)
            ax.set_ylim(-0.02, max(1.0, y_max) if y_max > 0 else 1.0)

    def _desenhar_barras(self, ax, df, col_metrica):
        is_en = getattr(self, "var_lang_en", None) and self.var_lang_en.get()
        lbl_better = " (lower = better)" if is_en else " (menor = melhor)"
        lbl_title = "Comparison by Compressor: " if is_en else "Comparação por Compressor: "

        ax.set_xlabel("Compressor", fontsize=10)
        if col_metrica == "cts":
            ax.set_ylabel("CTS (Custom Topological Score)" + lbl_better, fontsize=10)
        else:
            ax.set_ylabel(self.cb_metrica.get() + lbl_better, fontsize=10)
        ax.set_title(lbl_title + self.cb_metrica.get(), fontsize=11, fontweight="bold")
        ax.grid(True, axis="y", linestyle=":", alpha=0.6)

        if df.empty:
            ax.text(0.5, 0.5, "Sem dados para exibir", ha="center", va="center", color="gray")
            return

        ordem_compressores = ["bzip2", "zstd", "brotli", "ppmd", "Ficd/21g", "IRR-C/aFicd"]
        compressores_presentes = [c for c in ordem_compressores if c in df["compressor"].values]
        if not compressores_presentes:
            compressores_presentes = sorted(df["compressor"].unique())

        algoritmos_ordem = ["Neighbor-Joining", "UPGMA", "Complete Linkage"]
        cores_algo = [COR_NJ, COR_UPGMA, COR_CL]

        largura_barra = 0.25
        x_pos = []
        x_labels = []
        for i, comp in enumerate(compressores_presentes):
            base = i * (len(algoritmos_ordem) + 0.5)
            x_pos.append(base + largura_barra)
            x_labels.append(comp)

        import numpy as np
        for j, algo in enumerate(algoritmos_ordem):
            offset = j * largura_barra
            valores = []
            for comp in compressores_presentes:
                sub = df[(df["compressor"] == comp) & (df["algoritmo"] == algo)]
                if not sub.empty and sub[col_metrica].notna().any():
                    valores.append(sub[col_metrica].mean())
                else:
                    valores.append(np.nan)

            x_pos_algo = [x_pos[i] + offset - largura_barra for i in range(len(compressores_presentes))]
            ax.bar(x_pos_algo, valores, width=largura_barra, label=algo,
                   color=cores_algo[j], edgecolor="black", linewidth=0.5)

        valores_validos = df[col_metrica].dropna()
        if not valores_validos.empty:
            melhor = valores_validos.min()
            lbl_best = "Best" if (getattr(self, "var_lang_en", None) and self.var_lang_en.get()) else "Melhor"
            ax.axhline(melhor, color="#e74c3c", linestyle="--", linewidth=1.5,
                       label=f"{lbl_best} {self.cb_metrica.get()} = {melhor:.4f}")

        ax.set_xticks([pos for pos in x_pos])
        ax.set_xticklabels(x_labels, rotation=15, ha="right", fontsize=8)

        ax.legend(loc="best", fontsize=8, framealpha=0.8)
        ax.set_ylim(bottom=0)

    @staticmethod
    def _abreviatura_compressor(comp: str) -> str:
        if not comp:
            return "?"
        comp = comp.lower()
        if "bzip2" in comp:
            return "bzip2"
        if "zstd" in comp:
            return "zstd"
        if "brotli" in comp:
            return "brotli"
        if "ppmd" in comp:
            return "ppmd"
        if "ficd" in comp and "irr" in comp:
            return "aFicd"
        if "ficd" in comp:
            return "Ficd"
        if "replica" in comp:
            return "repl"
        if "asym" in comp:
            return "asym"
        if "21g" in comp:
            return "21g"
        if "modelo" in comp:
            return "Mod"
        if "artigo" in comp:
            return "Art"
        return comp[:6]

    def _exportar_graficos(self):
        if self.df_atual is None or self.df_atual.empty:
            messagebox.showinfo("Sem dados", "Carregue os dados primeiro.")
            return

        destino = filedialog.asksaveasfilename(
            title="Exportar gráficos (serão salvos _scatter e _bar em .png e .svg)",
            defaultextension=".png",
            filetypes=[("Imagem PNG", "*.png")],
            initialfile="graficos_comparativos",
            initialdir=_initialdir(PASTA_PADRAO_ARVORES),
        )
        if not destino:
            return

        base = os.path.splitext(destino)[0]
        try:
            va.exportar_figura(self.fig_scatter, base + "_scatter")
            va.exportar_figura(self.fig_bar, base + "_bar")
            messagebox.showinfo("Exportado", f"Gráficos salvos em:\n{base}_scatter.png/.svg\n{base}_bar.png/.svg")
        except Exception as e:
            messagebox.showerror("Erro ao exportar", str(e))

    def _exportar_tabela(self):
        if self.df_atual is None or self.df_atual.empty:
            messagebox.showinfo("Sem dados", "Carregue os dados primeiro.")
            return

        destino = filedialog.asksaveasfilename(
            title="Exportar tabela (Excel)",
            defaultextension=".xlsx",
            filetypes=[("Planilha Excel", "*.xlsx")],
            initialfile="dados_graficos",
            initialdir=_initialdir(PASTA_PADRAO_ARVORES),
        )
        if not destino:
            return

        try:
            self.df_atual.to_excel(destino, index=False)
            messagebox.showinfo("Exportado", f"Tabela salva em:\n{destino}")
        except ImportError:
            messagebox.showerror("Dependência ausente", "Instale openpyxl para exportar .xlsx")
        except Exception as e:
            messagebox.showerror("Erro ao exportar", str(e))


# ========================================================================= #
# Janela principal
# ========================================================================= #
class AppPrincipal(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Analisador Filogenético NCD — PIBIC (v5)")
        self.geometry("1320x900")

        notebook = ttk.Notebook(self)
        notebook.pack(fill="both", expand=True)

        self.aba_visualizador = AbaVisualizador(notebook)
        self.aba_metricas = AbaMetricas(notebook, callback_abrir_no_visualizador=self._abrir_no_visualizador)
        self.aba_comparador = AbaComparador(notebook)
        self.aba_graficos = AbaGraficos(notebook)

        notebook.add(self.aba_visualizador, text="  Visualizador de Árvores  ")
        notebook.add(self.aba_metricas, text="  Analisador de Métricas  ")
        notebook.add(self.aba_comparador, text="  Comparador de Árvores  ")
        notebook.add(self.aba_graficos, text="  📊 Gráficos Comparativos  ")

        self.notebook = notebook

        status = ttk.Frame(self)
        status.pack(side="bottom", fill="x")
        ttk.Separator(status, orient="horizontal").pack(fill="x")
        texto_status = f"Árvores: {PASTA_PADRAO_ARVORES}   |   Experimentais: {PASTA_EXPERIMENTAIS}   |   Relatórios: {PASTA_PADRAO_RELATORIO}"
        if not os.path.isdir(PASTA_PADRAO_ARVORES):
            texto_status = "⚠ Pasta padrão de árvores não encontrada — " + texto_status
        ttk.Label(status, text=texto_status, font=("Segoe UI", 8), foreground="#555555").pack(side="left", padx=6, pady=2)

    def _abrir_no_visualizador(self, nome_arquivo_nwk: str):
        caminho_arvores = os.path.join(PASTA_PADRAO_ARVORES, nome_arquivo_nwk)
        caminho_exp = os.path.join(PASTA_EXPERIMENTAIS, nome_arquivo_nwk)

        if os.path.exists(caminho_arvores):
            caminho = caminho_arvores
        elif os.path.exists(caminho_exp):
            caminho = caminho_exp
        else:
            caminho = filedialog.askopenfilename(
                title=f"Localize {nome_arquivo_nwk}",
                initialdir=_initialdir(PASTA_PADRAO_ARVORES),
                initialfile=nome_arquivo_nwk,
                filetypes=[("Newick tree", "*.nwk"), ("Todos os arquivos", "*.*")],
            )
            if not caminho:
                return

        self.aba_visualizador.carregar(caminho)
        self.notebook.select(self.aba_visualizador)



def main():
    app = AppPrincipal()
    app.mainloop()


if __name__ == "__main__":
    main()
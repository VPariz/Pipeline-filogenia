#!/usr/bin/env python3
"""
Script principal do Pipeline NCD.
Gera as matrizes de distancia usando os compressores (bzip2, zstd, brotli, ppmd).

Tem 2 modos aqui dentro:
1) replica_bat: Modo antigo so pra testar e comparar com os resultados velhos
   imitando o jeito que os scripts originais formatavam o texto (com as
   quebras de linha \r\n etc), mas essa parte nao entra no artigo final
2) linear_puro: A forma certa (e que gerou os resultados do BSB). Le a sequencia 
   inteira limpinha em uppercase, junta com um '$' e guarda o valor NCD real em float.
"""

import sys
import time
import shutil
import subprocess
import tempfile
import traceback
from pathlib import Path
from itertools import product as iter_product

import pandas as pd
from tqdm import tqdm

# Caminho das sequencias genomicas publicadas
_THIS_DIR = Path(__file__).parent.resolve()
INPUT_DIR  = (_THIS_DIR / ".." / ".." / "dados" / "genomas").resolve()
OUTPUT_DIR = Path("resultados_ncd")
VALID_EXTS = {".fasta", ".fa", ".fna", ".txt"}

REPLICA_SEP  = b"$\r\n"
OFFICIAL_SEP = b"$"


# =================== COMPRESSORES ===================

_BZIP2_EXE = shutil.which("bzip2")

def bz2_size(data: bytes) -> int:
    """bzip2 externo via subprocess; fallback para módulo bz2 do Python."""
    if _BZIP2_EXE:
        tmp_path = None
        bz2_path = None
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".tmp") as tmp:
                tmp_path = Path(tmp.name)
                tmp.write(data)
            bz2_path = tmp_path.with_suffix(".tmp.bz2")
            subprocess.run(
                [_BZIP2_EXE, "-fk", str(tmp_path)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return bz2_path.stat().st_size
        finally:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
            if bz2_path and bz2_path.exists():
                bz2_path.unlink()
    import bz2
    return len(bz2.compress(data, compresslevel=9))

try:
    import zstandard as zstd
    _ZSTD = zstd.ZstdCompressor(level=22)
    def zstd_size(data: bytes) -> int: return len(_ZSTD.compress(data))
    ZSTD_OK = True
except ImportError:
    ZSTD_OK = False

try:
    import brotli
    def brotli_size(data: bytes) -> int: return len(brotli.compress(data, quality=11))
    BROTLI_OK = True
except ImportError:
    BROTLI_OK = False

try:
    import pyppmd
    def ppmd_size(data: bytes) -> int:
        for kw in [{"max_order": 16, "mem_size": 64 << 20},
                   {"max_order": 8,  "mem_size": 16 << 20}, {}]:
            try:
                return len(pyppmd.compress(data, **kw))
            except TypeError:
                continue
        raise RuntimeError("pyppmd: nenhuma assinatura compatível.")
    PPMD_OK = True
except ImportError:
    PPMD_OK = False

# Todos os compressores disponíveis (mesma lista para ambas as abordagens)
COMPRESSORS: dict = {"bzip2": bz2_size}
if ZSTD_OK:   COMPRESSORS["zstd"]   = zstd_size
if BROTLI_OK: COMPRESSORS["brotli"] = brotli_size
if PPMD_OK:   COMPRESSORS["ppmd"]   = ppmd_size


# =================== LEITURA ===================

def read_replica(filepath: Path) -> bytes:
    """
    Simula FOR /F + echo %%L do Windows CMD:
      - Lê bytes brutos; normaliza terminadores para \\n
      - Descarta linhas vazias (FOR /F) e cabeçalhos '>'
      - Força \\r\\n ao final de cada linha sobrevivente (echo %%L)
    Resultado: linha1\\r\\nlinha2\\r\\n...últimaLinha\\r\\n
    """
    raw = filepath.read_bytes()
    lines = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n").split(b"\n")
    filtered = [ln for ln in lines if ln and not ln.startswith(b">")]
    if not filtered:
        raise ValueError(f"'{filepath.name}': sem sequência após filtro.")
    return b"\r\n".join(filtered) + b"\r\n"


def read_linear(filepath: Path) -> bytes:
    """
    DNA linear contínuo: remove cabeçalhos, quebras de linha e espaços,
    converte para maiúsculas. Sem filtros IUPAC.
    """
    parts = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(">") or s.startswith(";"):
                continue
            parts.append(s.replace(" ", "").upper())
    seq = "".join(parts)
    if not seq:
        raise ValueError(f"'{filepath.name}': sem sequência após limpeza.")
    return seq.encode("ascii")


# =================== NCD ===================

def ncd_integer(cx: int, cy: int, cxy: int) -> float:
    """Truncamento inteiro 32-bit do CMD: (10000*(cxy-min))//max / 10000"""
    mn, mx = min(cx, cy), max(cx, cy)
    if mx == 0:
        return 0.0
    return ((10000 * (cxy - mn)) // mx) / 10000.0


def ncd_real(cx: int, cy: int, cxy: int) -> float:
    """NCD float real de alta precisão."""
    mx = max(cx, cy)
    if mx == 0:
        return 0.0
    return (cxy - min(cx, cy)) / mx


# =================== MATRIZ ===================

def compute_matrix(species, sequences, comp_fn, cache,
                   sep: bytes, integer_ncd: bool, skip_diagonal: bool) -> list:
    n = len(species)
    mat = [[float("nan")] * n for _ in range(n)]
    seqs = [sequences[sp] for sp in species]

    for i, j in iter_product(range(n), range(n)):
        if skip_diagonal and i == j:
            continue
        cx, cy = cache[species[i]], cache[species[j]]
        if cx is None or cy is None:
            continue
        try:
            cxy = comp_fn(seqs[i] + sep + seqs[j])
            mat[i][j] = ncd_integer(cx, cy, cxy) if integer_ncd else ncd_real(cx, cy, cxy)
        except Exception as exc:
            tqdm.write(f"    [ERRO] NCD({species[i]},{species[j]}): {exc}")
    return mat


# =================== SALVAMENTO ===================

def save_matrix(mat, species, path: Path, decimals: int, skip_diagonal: bool):
    fmt = f".{decimals}f"
    with open(path, "w", encoding="utf-8") as f:
        for i, sp_i in enumerate(species):
            for j, sp_j in enumerate(species):
                if skip_diagonal and i == j:
                    continue
                val = mat[i][j]
                f.write(
                    f"NCD({sp_i},{sp_j})="
                    + ("nan" if pd.isna(val) else format(val, fmt))
                    + "\n"
                )


# =================== RUNNER ===================

def run_approach(label, suffix, species, sequences,
                 sep, integer_ncd, skip_diagonal, decimals):
    """Executa todos os compressores para uma abordagem e salva os resultados."""
    sep_repr = repr(sep)[2:-1]  # exibe b"$\r\n" → $\r\n
    diag_info = "sem diagonal" if skip_diagonal else "diagonal calculada"
    print(f"\n{'=' * 68}")
    print(f"  {label}")
    print(f"  sep=b\"{sep_repr}\"  |  {'inteiro 4d' if integer_ncd else 'float 6d'}  |  {diag_info}")
    print(f"{'=' * 68}")
    print(f"  {len(species)} espécies × {len(COMPRESSORS)} compressor(es): "
          f"{', '.join(COMPRESSORS)}")

    for comp_name, comp_fn in COMPRESSORS.items():
        print(f"\n  -- {comp_name.upper()} --")

        # Cache C(x)
        cache = {}
        for sp in tqdm(species, desc=f"  cache {comp_name}", unit="sp", leave=False):
            try:
                cache[sp] = comp_fn(sequences[sp])
                tqdm.write(f"    C({sp}) = {cache[sp]:,} bytes")
            except Exception as exc:
                tqdm.write(f"    [ERRO cache] {sp}: {exc}")
                cache[sp] = None

        mat = compute_matrix(species, sequences, comp_fn, cache,
                             sep, integer_ncd, skip_diagonal)

        out_path = OUTPUT_DIR / f"ncd-{comp_name}-{suffix}.txt"
        save_matrix(mat, species, out_path, decimals, skip_diagonal)
        tqdm.write(f"  -> {out_path.name}")


# =================== MAIN ===================

def main():
    t0 = time.time()
    print("=" * 68)
    print("  NCD Pipeline v13  —  pasta única  |  2 abordagens  |  todos os compressores")
    print("=" * 68)

    if not INPUT_DIR.exists():
        print(f"\n[ERRO] Pasta '{INPUT_DIR}' não encontrada.")
        sys.exit(1)

    fasta_files = sorted(
        f for f in INPUT_DIR.iterdir()
        if f.suffix.lower() in VALID_EXTS and f.is_file()
    )
    if not fasta_files:
        print(f"\n[ERRO] Nenhum arquivo FASTA em '{INPUT_DIR}'.")
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n  Entrada  : {INPUT_DIR.resolve()}")
    print(f"  Saída    : {OUTPUT_DIR.resolve()}")
    print(f"  Arquivos : {len(fasta_files)}")
    for f in fasta_files:
        print(f"    • {f.name}")
    print(f"\n  bzip2 externo : "
          f"{'SIM — ' + str(_BZIP2_EXE) if _BZIP2_EXE else 'NÃO (fallback Python bz2)'}")
    print(f"  Compressores  : {', '.join(COMPRESSORS)}")

    # ── Leitura ────────────────────────────────────────────────────────
    seqs_replica, seqs_linear = {}, {}

    print("\n  [Leitura]")
    for fpath in tqdm(fasta_files, desc="  réplica", unit="arq"):
        try:
            seqs_replica[fpath.name] = read_replica(fpath)
        except Exception as exc:
            tqdm.write(f"  ✗ réplica {fpath.name}: {exc}")

    for fpath in tqdm(fasta_files, desc="  linear ", unit="arq"):
        try:
            seqs_linear[fpath.name] = read_linear(fpath)
        except Exception as exc:
            tqdm.write(f"  ✗ linear  {fpath.name}: {exc}")

    species_replica = sorted(seqs_replica)
    species_linear  = sorted(seqs_linear)

    # ── Abordagem 1: replica_bat ───────────────────────────────────────
    if seqs_replica:
        run_approach(
            label        = "ABORDAGEM 1 — replica_bat",
            suffix       = "replica_bat",
            species      = species_replica,
            sequences    = seqs_replica,
            sep          = REPLICA_SEP,   # b"$\r\n"  (type A já termina em \r\n)
            integer_ncd  = True,
            skip_diagonal= True,
            decimals     = 4,
        )

    # ── Abordagem 2: linear_puro ───────────────────────────────────────
    if seqs_linear:
        run_approach(
            label        = "ABORDAGEM 2 — linear_puro",
            suffix       = "linear_puro",
            species      = species_linear,
            sequences    = seqs_linear,
            sep          = OFFICIAL_SEP,  # b"$"
            integer_ncd  = False,
            skip_diagonal= False,
            decimals     = 6,
        )

    # ── Resumo ─────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    arquivos_gerados = sorted(OUTPUT_DIR.glob("ncd-*.txt"))

    print(f"\n{'=' * 68}")
    print("  RESUMO")
    print(f"{'=' * 68}")
    print(f"  Espécies       : {len(fasta_files)}")
    print(f"  Compressores   : {', '.join(COMPRESSORS)}")
    print(f"  Arquivos gerados ({len(arquivos_gerados)}):")
    for p in arquivos_gerados:
        print(f"    • {p.name}")
    print(f"  Tempo total    : {elapsed:.2f} s")
    print(f"{'=' * 68}")
    print("  Pipeline v13 concluído.")
    print(f"{'=' * 68}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[INFO] Interrompido.")
        sys.exit(0)
    except Exception as exc:
        print(f"\n[ERRO FATAL] {exc}")
        traceback.print_exc()
        sys.exit(1)

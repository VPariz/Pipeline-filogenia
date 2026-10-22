# Pairwise NCD Distance Matrices

This directory contains the canonical pairwise distance matrices computed via Normalized Compression Distance (NCD) for the 20 mammalian mitochondrial genomes benchmark dataset.

## File Format

Each matrix file is stored in plain text (`.txt`), with each line representing the pairwise distance between two genomes:

```text
NCD(genome_A.txt, genome_B.txt) = float_value
```

The metric is symmetric with $NCD(x, y) \in [0, 1 + \varepsilon]$, computed according to the formal definition by Cilibrasi & Vitányi (2005):

$$NCD(x, y) = \frac{C(xy) - \min(C(x), C(y))}{\max(C(x), C(y))}$$

where $C(s)$ denotes the compressed size in bytes of string $s$, and $xy$ represents the concatenation of the two genome sequences with a standard `$` separator.

---

## Included Matrices

| File | Compressor / Algorithm | Configuration / Parameters | Taxa |
|---|---|---|:---:|
| `zstd.txt` | **Zstandard** | Compression level `22` (maximum level in `zstandard`) | 20 |
| `brotli.txt` | **Brotli** | Quality level `11` (maximum quality level in `brotli`) | 20 |
| `bzip2.txt` | **bzip2** | Level `-9` (900 KB block size, pure continuous sequence) | 20 |
| `ppmd.txt` | **PPMD (PPMII)** | Model order `16`, memory limit `64 MB` (`pyppmd`) | 20 |
| `aficd.txt` | **IRR-C** | Admissible CFG heuristic selecting rules with highest compression | 21* |

*\* Note: `aficd.txt` includes the 20 mammalian taxa plus an artificial random control outgroup (`randomFasta-N16k.txt`).*

---

## Loading in Python

```python
import re
import numpy as np
import pandas as pd

def load_ncd_matrix(filepath):
    distances = {}
    taxa = set()
    pattern = re.compile(r"NCD\(([^,]+),\s*([^)]+)\)\s*=\s*([0-9.]+)")
    
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                s1, s2, dist = m.group(1).strip(), m.group(2).strip(), float(m.group(3))
                distances[(s1, s2)] = dist
                taxa.update([s1, s2])
                
    taxa = sorted(taxa)
    n = len(taxa)
    mat = np.zeros((n, n))
    for i, t1 in enumerate(taxa):
        for j, t2 in enumerate(taxa):
            mat[i, j] = distances.get((t1, t2), distances.get((t2, t1), 0.0))
            
    return pd.DataFrame(mat, index=taxa, columns=taxa)

df_zstd = load_ncd_matrix("zstd.txt")
print(df_zstd.shape)
```

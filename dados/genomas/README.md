# Mitochondrial Genomes Benchmark Dataset (mtDNA)

This directory contains the complete mitochondrial genomes (mtDNA) of 20 mammalian species representing major evolutionary lineages, originally compiled by Li et al. (2001, 2004) and Otu & Sayood (2003).

All sequences were retrieved from the **NCBI GenBank** repository (Sayers et al., 2022).

---

## Taxa and Accession Numbers

| File Name | Common Name (EN) | Nome Comum (PT) | Scientific Name | Clade / Taxonomic Order | GenBank Accession | Length (bp) |
|---|---|---|---|---|---|:---:|
| `blue_whale.txt` | Blue Whale | Baleia-azul | *Balaenoptera musculus* | Cetartiodactyla | [`NC_001601.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001601.1) | 16,402 |
| `bornean_orangutan.txt` | Bornean Orangutan | Orangotango-de-bornéu | *Pongo pygmaeus* | Primates | [`NC_001646.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001646.1) | 16,389 |
| `cat.txt` | Domestic Cat | Gato | *Felis catus* | Carnivora | [`NC_001700.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001700.1) | 17,009 |
| `chimpanzee.txt` | Chimpanzee | Chimpanzé | *Pan troglodytes* | Primates | [`NC_001643.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001643.1) | 16,554 |
| `common_wallaroo.txt` | Common Wallaroo | Canguru-euro | *Macropus robustus* | Diprotodontia (Marsupial) | [`NC_001794.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001794.1) | 16,896 |
| `cow.txt` | Cow | Vaca | *Bos taurus* | Cetartiodactyla | [`NC_006853.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_006853.1) | 16,338 |
| `finback_whale.txt` | Finback Whale | Baleia-comum | *Balaenoptera physalus* | Cetartiodactyla | [`NC_001321.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001321.1) | 16,398 |
| `gibbon.txt` | Lar Gibbon | Gibão | *Hylobates lar* | Primates | [`NC_002082.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_002082.1) | 16,472 |
| `gorilla.txt` | Western Gorilla | Gorila | *Gorilla gorilla* | Primates | [`NC_001645.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001645.1) | 16,364 |
| `grey_seal.txt` | Grey Seal | Foca-cinzenta | *Halichoerus grypus* | Carnivora | [`X72004.1`](https://www.ncbi.nlm.nih.gov/nuccore/X72004.1) | 16,797 |
| `harbor_seal.txt` | Harbor Seal | Foca-comum | *Phoca vitulina* | Carnivora | [`NC_001325.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001325.1) | 16,826 |
| `horse.txt` | Horse | Cavalo | *Equus caballus* | Perissodactyla | [`AB859014.1`](https://www.ncbi.nlm.nih.gov/nuccore/AB859014.1) | 16,616 |
| `house_mouse.txt` | House Mouse | Camundongo | *Mus musculus* | Rodentia | [`D83491.1`](https://www.ncbi.nlm.nih.gov/nuccore/D83491.1) | 16,300 |
| `human.txt` | Human | Humano | *Homo sapiens* | Primates | [`AP009475.1`](https://www.ncbi.nlm.nih.gov/nuccore/AP009475.1) | 16,568 |
| `north_american_opossum.txt` | North American Opossum | Gambá-da-virgínia | *Didelphis virginiana* | Didelphimorphia (Marsupial) | [`NC_001610.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001610.1) | 17,084 |
| `norway_rat.txt` | Norway Rat | Rato-marrom | *Rattus norvegicus* | Rodentia | [`MZ014548.1`](https://www.ncbi.nlm.nih.gov/nuccore/MZ014548.1) | 16,310 |
| `platypus.txt` | Platypus | Ornitorrinco | *Ornithorhynchus anatinus* | Monotremata | [`NC_000891.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_000891.1) | 17,019 |
| `pygmy_chimpanzee.txt` | Pygmy Chimpanzee (Bonobo) | Bonobo | *Pan paniscus* | Primates | [`NC_001644.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001644.1) | 16,563 |
| `sumatran_orangutan.txt` | Sumatran Orangutan | Orangotango-de-sumatra | *Pongo abelii* | Primates | [`NC_002083.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_002083.1) | 16,499 |
| `white_rhino.txt` | White Rhinoceros | Rinoceronte-branco | *Ceratotherium simum* | Perissodactyla | [`NC_001808.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001808.1) | 16,832 |

---

## Sequence Format

Each `.txt` file contains the complete nucleotide sequence (IUPAC A, C, G, T) in uppercase, formatted as raw continuous text. When processed by the compression pipeline, sequences are read into memory and whitespace/linebreaks are stripped to ensure exact, unpadded compression calculation.

## References

1. **Li, M. et al. (2001)**. An information-theoretic metric and its application to phylogeny. *Journal of the ACM*, 51(6), 878-894.
2. **Li, M. et al. (2004)**. The similarity metric. *IEEE Transactions on Information Theory*, 50(12), 3250-3264.
3. **Otu, H. H., & Sayood, K. (2003)**. A new sequence distance measure for phylogenetic tree construction. *Bioinformatics*, 19(16), 2122-2130.
4. **Sayers, E. W. et al. (2022)**. Database resources of the National Center for Biotechnology Information. *Nucleic Acids Research*, 50(D1), D20-D26.

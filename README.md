# TidyTool

**TidyTool** is a desktop application for analyzing, filtering, and cleaning biological sequences in FASTA format. It provides a powerful and intuitive **PyQt6-based GUI** for:

- Removing duplicate sequences,
- Filtering sequences by name, length, or similarity,
- Parsing and visualizing BLAST results,
- Displaying histograms for score, E-value, and sequence length.

---

## Features

- Load multiple FASTA files at once
- Remove duplicates (by name, sequence, or both)
- Filter sequences based on:
  - Name (partial or full match),
  - Sequence similarity (custom threshold),
  - Length range (min/max),
  - Start codon (`M` for proteins, `ATG` for genes),
  - Score and E-value (parsed from BLAST headers)
- Interactive histogram display:
  - Lengths
  - Score
  - E-values — switchable by button
- Export selected sequences to a FASTA file
- Right-click on a sequence → see full FASTA in a copyable popup

---

## Requirements

- Python `>=3.7`
- `PyQt6`
- `matplotlib`
- `sqlite3` (built into Python)

Install dependencies:
```bash
pip install -r requirements.txt
```

---

## How to Run

1. Launch the app:
```bash
python main.py
```

2. In the **CLEANER** tab:
   - Click **Load File** and select `.fasta` files,
   - Use **Clean** to remove duplicates, or **Filter** to search,
   - Select sequences from the list and export/save them.

3. In the **BLAST** tab:
   - Choose a BLAST database,
   - Run BLAST search with selected sequences,
   - Open and visualize the results with interactive filters.

---

## Folder Structure

```
.
├── Libraries
├	├── blast_parser.py         # BLAST XML parsing
├	├── blast_utils.py          # BLAST execution
├	├── blast_view.py           # BLAST plot view
├	├── draw_utils.py           # plotting logic
├	├── fasta_utils.py          # FASTA logic
├	├── gui.py                  # main GUI
├	├── move_utils.py           # list management helpers
├	├── thread_utils.py         # background threads
├── BLAST
├	├── blastn.exe
├	├── blastp.exe
├	├── blastx.exe
├	├── makeblastdb.exe
├	├── ncbi-vdb-md.dll
├	├── nghttp2.dll
├── BAZA
├	├── database_here.txt
├── Resources
├	├── tidytool.png
├── main.py
├── requirements.txt
└── README.md
```

---

## Sample Data

You can use any FASTA-formatted files containing protein or gene sequences.
We recommend UniProt, Ensembl, or NCBI RefSeq for testing.

---

## Authors

- Aleksandra Liszka
- Artur Stołowski

---

## License

This project is open-source for educational and research purposes.  
There is no guarantee of biological correctness — use at your own risk 😉

# TidyTool

**TidyTool** is a desktop application for analyzing, filtering, and cleaning biological sequences in FASTA format. It provides a powerful and intuitive **PyQt6-based GUI** for:

- Removing duplicate sequences,
- Filtering sequences by name, length, or similarity,
- Parsing and visualizing BLAST results,
- Displaying histograms for score, E-value, and sequence length.

---
![ScreenShot](resources/b_screen.png)
![ScreenShot2](resources/c_screen.png)
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

Windows:
Install dependencies:
```bash
pip install -r requirements.txt
```
Linux:
```bash
sudo apt install python3-pyqt6 python3-matplotlib
```
Mac:
```bash
pip3 install -r requirements.txt
```
---

## How to Run

1. Launch the app:

Windows:
```bash
python main.py
```
Linux:
```bash
python3 main.py
```
Mac:
```bash
python3 main.py
```

2. In the **CLEANER** tab:
   - Click **Load File** and select `.fasta` files,
   - Use **Clean** to remove duplicates, or **Filter** to search,
   - Select sequences from the list and export/save them.

3. In the **BLAST** tab:
   - Choose a BLAST database,
   - Run BLAST search with selected sequences,
   - Open and visualize the results.

---

## Folder Structure

```
.
├── libraries
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
├	├── blastn
├	├── blastp.exe
├	├── blastp
├	├── blastx.exe
├	├── blastx
├	├── tblastx.exe
├	├── tblastx
├	├── tblastn.exe
├	├── tblastn
├	├── makeblastdb.exe
├	├── makeblastdb
├	├── blast_formatter.exe
├	├── blast_formatter
├	├── ncbi-vdb-md.dll
├	├── nghttp2.dll
├── database
├	├── database_here.txt
├── resources
├	├── tidytool.png
├── main.py
├── requirements.txt
├── LICENCE
├── NCBI_LICENCE
└── README.md
```

---

## Sample Data

You can use any FASTA-formatted files containing protein or gene sequences.
We recommend UniProt, Ensembl, or NCBI RefSeq for testing.

---

## Authors

- Aleksandra Liszka
- Aleksandra Marcisz
- Artur Stołowski

---

This application uses BLAST+ tools developed by the National Center for Biotechnology Information (NCBI).
https://blast.ncbi.nlm.nih.gov/

## License

MIT License

Copyright (c) 2025 Aleksandra Liszka, Artur Stołowski, Aleksandra Marcisz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

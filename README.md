# BLASTnBRUSH

**BLASTnBRUSH** is a user-friendly, graphical application designed to facilitate working with sequence databases, especially for researchers without advanced bioinformatics skills. It provides a powerful and intuitive **PyQt6-based GUI** for:

- Removing duplicate sequences,
- Filtering sequences by name, length, or similarity,
- Parsing and visualizing BLAST results,
- Displaying histograms for score, E-value, and sequence length.

---
![ScreenShot](resources/b_screen.png)
![ScreenShot2](resources/c_screen.png)
---

## Features
**BLAST+**:
- Load multiple FASTA files at once for creating BLAST database
- Run BLASTP, BLASTN, BLASTX, TBLASTN or TBLASTX
- Show Alignments (BLAST output file, -outfrm 0)
- Show Alignment coverage plots
- Export FASTA file with BLAST headers for further analyse in "Clean" tab

**Cleaner**:
- Load multiple FASTA files at once
- Remove duplicates (by name, sequence, or both)
- Filter sequences based on:
  - Name (partial or full match),
  - Sequence similarity (custom threshold),
  - Length range (min/max),
  - Start codon (`M` for proteins, `ATG` for genes),
  - Score, E-value, Alignment Length, Identity, Similarity (parsed from BLAST headers)
- Interactive histogram display:
  - Length
  - Score
  - E-value
  - Alignment Length
  - Identity
  - Similarity
- Export selected sequences to a FASTA file
- Right-click on a sequence → see full FASTA in a copyable popup

---

## Requirements

- Python `>=3.12`
- `PyQt6`
- `matplotlib`
- `sqlite3` (built into Python)


### **Windows:**
**- Python version:**
  - Download Python installer:
    https://www.python.org/downloads/windows/
    Choose Python 3.12.x and download 64-bit version.

  - Run installer
    SELECT "Add Python 3.12 to PATH" (very important).
  - Click "Install Now"
  - Download TidyTool.zip from Releases
  - Extract to desired folder
  - Open TidyTool folder
  - Press the right mouse button and select "Open in Terminal"
  - Paste commands:
  ```bash
  pip install -r requirements.txt
  ```
  - After installation type:
  ```bash
  python main.py
  ```
**- Exe version:**
  -  Download TidyTOOL Windows.zip from Releases
  -  Extract do desired folder
  -  Open TidyTool folder
  -  Run TidyTool.exe

**Important Notice:**
This application might be falsely detected as a virus by some antivirus software. This is a false positive, which is a common issue for apps generated from Python using tools like auto-py-to-exe or PyInstaller.


### **Linux:**
- Download TidyTool.zip from Releases
- Extract to desired folder
- Open TidyTool folder
- Right-click inside folder → "Open in Terminal"
- Paste command:
```bash
sudo apt install python3-pyqt6 python3-matplotlib
```
- Run TidyTOOL by command:
```bash
python3 main.py
```


### **Mac:**
- Download Python installer:
https://www.python.org/downloads/macos/
Choose Python 3.12.x and download 64-bit version:
"macOS 64-bit universal2 installer"
- Run installer
- Click "Install Now"
- Download TidyTool.zip from Releases
- Extract to desired folder
- Open TidyTool folder
- Run Terminal (Command (⌘) + Space)
- Navigate to downloaded TidyTool folder
(eg./Users/username/Downloads/TidyTool)
- Paste command:
```bash
pip3 install -r requirements.txt
```
- Run application:
```bash
python3 main.py
```
---


## How to Use

1. In the **BLAST** tab:
   - Choose a FASTA file to create BLAST database (Protein or Nucleotide),
   - Paste query sequences in Query box,
   - Run suitable BLAST search (BLASTP, BLASTN, BLASTX, TBLAXTN, TBLASTN),
   - Click on **Show Alignment** to open full BLAST result,
   - Select query in Query List and click on **View** to open Alignment coverage plot,
   - Click on **Save sequences** to create FASTA file (contains headers combined with BLAST data and belonging sequences).
2. In the **Cleaner** tab:
   - Click **Choose FASTA File** and select `.fasta` files,
   - You can preview the files you selected by clicking the **Show Files** button,
   - Select the method for removing duplicates:
     - Names (removes sequences with the same headers),
     - Sequences (removes identical sequences),
     - None selected (does not remove duplicates),
   - Click **Clean/Analyse** to start analyzing the selected files,
   - The headers will appear in the **Genes** window,
   if a file containing integrated BLAST information is loaded, the Score and E-value histograms will be generated automatically.
   - You can choose which histogram to display in each of the three panels by clicking the corresponding buttons,
   - You can choose which filter to apply by clicking the corresponding checkbox next to its name,
   - After entering the filter parameters, click the **Filter** button to display only the sequences that match the selected filters,
   - Selected genes can be highlighted in the **Genes** list and moved to the **Selected names** list using the **ADD** button,
   - By clicking the **Save FASTA** button, you save the genes from the **Selected names** list as a new FASTA file.


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

All databases used in testing the program can be download from following sources.

**CAZy Database** - Carbohydrate Active Enzymes database, http://www.cazy.org/

Cantarel BL, Coutinho PM, Rancurel C, Bernard T, Lombard V, Henrissat B (2009) The Carbohydrate-Active EnZymes
database (CAZy) : an expert resource for Glycogenomics. Nucleic Acids Res 37:D233-238 [PMID : 18838391].
```
https://bcb.unl.edu/dbCAN2/download/CAZyDB.07142024.fa
```

**TreeGenes** - https://treegenesdb.org/

Wegrzyn J.L., Staton M.A., Street N. R., Main D., Grau E., Herndon N., Buehler S., Falk T., Zaman S., Ramnath R., Richter
P., Sun L., Condon B., Almsaeed A., Chen M.,Mannapperuma C., Jung S., Ficklin S. Cyberinfrastructure to Improve Forest
Health and Productivity: The Role of Tree Databases in Connecting Genomes, Phenomes, and the Environment,
TreeGenes. Database, Volume 2019. doi:10.3389/fpls.2019.00813

Falk T., Herndon N., Grau E., Buehler S., Richter P., Zaman S., Baker E.M., Ramnath R., Ficklin S., Staton M., Feltus F.A.,
Jung S., Main D., Wegrzyn J.L. (2018). Growing and cultivating the forest genomics database, TreeGenes. Database,
Volume 2018 doi:10.1093/database/bay084
```
Abies alba
https://treegenesdb.org/FTP/Genomes/Abal/v1.1/annotation/Abal.1_1.pep.fa.gz
Picea abies
https://treegenesdb.org/FTP/Genomes/Paab/v1.0b/annotation/Paab.1_0b.pep.fa.gz
Populus tremula x tremuloides
https://treegenesdb.org/FTP/Genomes/Potl/v1.0/annotation/Potl.1_0.pep.fa.gz
Ginkgo Biloba
https://treegenesdb.org/FTP/Genomes/Gibi/v1.0/annotation/Gibi.1_0.pep.fa.gz
```

**OneKP** - https://db.cngb.org/onekp/

J.H. Leebens-Mack, M.S. Barker, E.J. Carpenter, M.K. Deyholos, M.A. Gitzendanner, S.W. Graham, I. Grosse, Z. Li, M.
Melkonian, S. Mirarab, ... G.K. Wong. Oct 2019. One thousand plant transcriptomes and the phylogenomics of green plants.
Nature 574: 679-685. PMID: 31645766. DOI: 10.1038/s41586-019-1693-2

E.J. Carpenter, N. Matasci, S. Ayyampalayam, S. Wu, J. Sun, J. Yu, F.R. Jimenez Vieira, C. Bowler, R.G. Dorrell, M.A.
Gitzendanner, ... G.K. Wong. Oct 2019. Access to RNA-sequencing data from 1173 plant species: The 1000 plant
transcriptomes initiative (1KP). GigaScience 8: giz126. PMID: 31644802. DOI: 10.1093/gigascience/giz126
```
https://ftp.cngb.org/pub/SciRAID/onekp/assemblies/DZQM-Pinus_radiata/DZQM-translated-protein.fa.gz
```

**TAIR database** - The Arabidopsis Information Resource, www.arabidopsis.org

Leonore Reiser, Erica Bakker, Sabarinath Subramaniam, Xingguo Chen, Swapnil Sawant, Kartik Khosa, Trilok Prithvi,
Tanya Z Berardini, The Arabidopsis Information Resource in 2024, Genetics, Volume 227, Issue 1, May 2024,
iyae027, https://doi.org/10.1093/genetics/iyae027
```
https://www.arabidopsis.org/download/file?path=Proteins/Araport11_protein_lists/Araport11_pep_20250411.gz
```

**Solgenomics** - https://solgenomics.sgn.cornell.edu/

Noe Fernandez-Pozo, Naama Menda, Jeremy D. Edwards, Surya Saha, Isaak Y. Tecle, Susan R. Strickler, Aureliano
Bombarely, Thomas Fisher-York, Anuradha Pujar, Hartmut Foerster, Aimin Yan, Lukas A. Mueller, The Sol Genomics
Network (SGN)—from genotype to phenotype to breeding, Nucleic Acids Research, Volume 43, Issue D1, 28 January 2015,
Pages D1036–D1041, https://doi.org/10.1093/nar/gku1195
```
https://solgenomics.net/ftp/genomes/Nicotiana_benthamiana/annotation/Niben101/Niben101_annotation.proteins.fasta.gz
```

---

## Authors

- Aleksandra Liszka
- Aleksandra Marcisz
- Artur Stołowski

---

This application uses BLAST+ tools developed by the National Center for Biotechnology Information (NCBI).
https://blast.ncbi.nlm.nih.gov/


# BLASTnBRUSH

**BLASTnBRUSH** is a user-friendly, graphical application designed to facilitate working with sequence databases, especially for researchers without advanced bioinformatics skills. It provides a powerful and intuitive **PyQt6-based GUI** for:

- Local BLAST analysis facilitating multiple database loading and custom-provided databases with multiple queries
- Visualizing BLAST results and alignment coverage for each queries
- Removing duplicate sequences
- Filtering both by the nucleotide/amino acid sequence and by sequence identifiers
- Filtering BLAST results (both user-performed BLAST and NCBI BLAST XML files) by sequence length, alignment length, E-value, bit score, identity, similarity, motif/domain search, and start codon detection
- Displaying histograms for BLAST metrics

---
---

## Features
**BLAST+**:
- Load multiple FASTA files simultaneously for use in BLAST searches
- Load BLAST XML result files and FASTA databases to link BLAST hits with the corresponding database sequences
- Add multiple queries at once 
- Run blastp, blastn, blastx, tblastn or tblastx
- Display alignments and BLAST results (BLAST output file, -outfrm 0)
- Show alignment coverage plots
- Export files with BLAST metrics and alignments
- Export FASTA files with BLAST metrics and sequences for further analysis in the "Cleaner" tab

**Cleaner**:
- Load multiple FASTA files to be processed
- Interactive histogram display for BLAST metrics (sequence length, Bit score, E-value, Alignment Length, Identity, Similarity)
- Remove duplicates from BLAST results (by name, sequence, or both)
- Filter and extract sequences based on:
  - Name (partial or full match),
  - domain/motifs similarity (partial or full match),
  - Sequence length (min/max),
  - Start codon presence (start metionin in proteins and ATG codon for genes),
  - BLAST metrics (Bit score, E-value, Alignment Length, Identity, Similarity (from XML file))
  - Right-click a sequence name in the list to view and copy its full FASTA format
- Export selected sequences to a FASTA file

---
# Instalation

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
  - Download BLASTnBRUSH.zip from Releases
  - Extract to desired folder
  - Open BLASTnBRUSH folder
  - Run "Install Python libraries (Windows).bat" - it will install PyQt6 and matplotlib libraries
  - Run "BLASTnBRUSH (Windows).bat" - it will open BLASTnBRUSH application in Python

**- Exe version:**
  -  Download BLASTnBRUSH Portable.zip from Releases
  -  Extract do desired folder
  -  Open BLASTnBRUSH folder
  -  Run BLASTnBRUSH.exe

**Important Notice:**
This application might be falsely flagged as a virus by some antivirus software. This is a false positive, which is a common issue for apps generated from Python using tools like auto-py-to-exe or PyInstaller.


### **Linux:**
- Download BLASTnBRUSH.zip from Releases
- Extract to desired folder
- Open BLASTnBRUSH folder
- Right-click on "Install Python libraries (Linux).sh" and select Run as Program
- Right-click on "BLASTnBRUSH (Linux).sh" and select Run as Program

### **Mac:**
- Download Python installer:
https://www.python.org/downloads/macos/
Choose Python 3.12.x and download 64-bit version:
"macOS 64-bit universal2 installer"
- Run installer
- Click "Install Now"
- Download BLASTnBRUSH.zip from Releases
- Extract to desired folder
- Open BLASTnBRUSH folder
- Press (Command (⌘) + Space), type Terminal in search window and run it
- Navigate to downloaded BLASTnBRUSH folder
(eg.cd Downloads/BLASTnBRUSH)
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
   - Choose a FASTA file(s) to create BLAST database(s) (Protein or Nucleotide),
   - Paste query sequence(s) in Query box,
   - Run suitable BLAST search (BLASTP, BLASTN, BLASTX, TBLAXTN, TBLASTN),
   - Click on **Show Alignment** to open full BLAST result, if you wish save if from the file,
   - Select query in Query List on the right and click on **View** to open Alignment coverage plot for each query. The number above aligement coverage informs about aligement lenght
   - Click on **Save sequences** to create FASTA file (contains seqience queries, combined with BLAST data and belonging amino acids /nucleotide sequences).
   - To browse and later filter BLAST results precomputed with NCBI: first, load the database with full sequences by clicking **Load database**, then import the BLAST XML by clicking **Import BLAST results XML**. After that, click **Save sequences** to bind the XML with full-length sequences, making them ready for screening with Cleaner.
2. In the **Cleaner** tab:
   - Click **Choose FASTA File** and select `.fasta` files as ready datasets ot file prepared in previous step,
   - You can preview the files you selected by clicking the **Show Files** button,
   - Select the method for removing duplicates:
     - By names (removes sequences with the same headers),
     - By sequences (removes identical sequences),
     - None selected (does not remove duplicates),
   - Click **Clean/Analyse** to start analyzing the selected files,
   - The headers will appear in the **Sequences** window,
 If a file containing integrated BLAST information is loaded, the user can select up to three histograms to display result distributions.
   - You can choose which histogram to display in each of the three panels by clicking the corresponding buttons,
   - You can choose which filter to apply by clicking the corresponding checkbox next to its name,
   - After entering the filter parameters, click the **Filter** button to display only the sequences that match the selected filters,
   - Selected sequences can be highlighted in the **sequences** list and moved to the **Selected names** list using the **ADD** button,
   - You can check for specific domains or motifs and their identity percentage when analyzing sequences,
   - By right-clicking on each sequence, you can view the full sequence for the selected name,
   - By clicking the Save FASTA button, you save the names and corresponding sequences from the Selected names list as a new FASTA file, ready to be analyzed in subsequent approaches.



# Sample Data

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


# Case Use

## 1.   **Performance and Scalability** 
   
   - The CAZy protein database was downloaded and loaded in the Cleaner tab by clicking **Load FASTA file**.
   - In the **Check duplicates section**, Sequences was selected, and the database was cleaned of duplicated sequences by clicking **Clean/Analyse**. (Since CAZy uses many sources of data repositories, the same sequence may be present under multiple names).
   - Protein names from the GT2 family belonging to *Arabidopsis thaliana* were copy-pasted into the Enter by name window, then Filter was clicked. GT2 family list is available on CAZy database under the link https://www.cazy.org/IMG/cazy_data/GT2.txt
      ![ScreenShot](resources/GT2_extract.png)
   - All sequences were selected by checking **Select All**, then clicked **Add**.
   - The database was saved by clicking **Save FASTA**, (database containing only the GT2 family from *Arabidopsis thaliana*).
   - In the BLAST tab, all extracted sequences were used as queries against databases for *Nicotiana benthamiana, Populus tremula x tremuloides, Abies alba, Pinus radiata, Picea abies, and Ginkgo biloba*
   - BLASTP was run
        ![ScreenShot](resources/BLASTP_results_GT2.png)
   - Results were saved by clicking **Save Sequences**
   - The saved results were loaded into the Cleaner tab. The BLAST results were deduplicated by **Clean/Analyse**.
    ![ScreenShot](resources/GT2_BLAST_cleaner.png)
   - Protein lengths were filtered between 100–800 amino acids, similarity thresholds set between 50%–100%, and alignment length between 100 and max. Then **Filter** was clicked.
   - **Select All** and **Add**  and **Save Fasta** were clicked
      ![ScreenShot](resources/GT2_BLAST_cleaned.png)
   - The results were ready to be exported to MEGA-X for phylogenetic analysis.

## 2. **Hypothesis-driven analysis**
   
   - The TAIR protein database was downloaded and loaded into the Cleaner section. **Clean/Analyse** was clicked without removing any duplicates.
   - In the **Enter by name field**, the text "CSLA, 1" was typed, and **Filter** was clicked to obtain all CSLA singular isoforms
      ![ScreenShot](resources/CSLA_extract.png)
   - All sequences were selected, **Add** clicked, and the **FASTA** file saved.
   - The obtained sequences were used as queries in the BLAST tab after loading databases for *Nicotiana benthamiana, Populus tremula x tremuloides, Abies alba, Pinus radiata, Picea abies, and Ginkgo biloba*.
   - BLASTP was run.
   - Results were saved by clicking **Save Sequences**.
      ![ScreenShot](resources/CSLA_BLAST.png)
   - The saved results were loaded into the Cleaner tab. The BLAST results were deduplicated. Protein lengths were filtered between 100–800 amino acids, similarity thresholds set between 50%–100%, and alignment length between 100 and max. Then **Filter** was clicked.
   - **Select all** and **ADD**
 ![ScreenShot](resources/CSLA_cleaned.png)
   - The results were ready to be expored to MEGA-X for phylogenetic analysis.
     

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
├	├── BLASTnBRUSH.png
├── BLASTnBRUSH (Linux).sh
├── BLASTnBRUSH (MacOS).command
├── BLASTnBRUSH (Windows).bat
├── Install Python libraries (Linux).sh
├── Install Python libraries (MacOS).command
├── Install Python libraries (Windows).bat
├── main.py
├── requirements.txt
├── LICENCE
├── NCBI_LICENCE
└── README.md
```

---



---

## Authors

- Aleksandra Liszka
- Aleksandra Marcisz
- Artur Stołowski

---

This application uses BLAST+ tools developed by the National Center for Biotechnology Information (NCBI).
https://blast.ncbi.nlm.nih.gov/


# This file is part of TidyTool
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the MIT License

from blast_parser import parse_blast_output
import subprocess
from PyQt6.QtWidgets import QFileDialog, QMessageBox




###################### CREATING DATABASE ########################################
#   LOAD PREVIOUS
def load_prev_database(widget, label_to_update):
    merged_path = "merged_database.fasta"
    try:
        import os
        
        label_to_update.setText("merged_database.fasta")
        db_name = os.path.splitext(os.path.basename(merged_path))[0]
        QMessageBox.information(widget, "Success", f"Database loaded: {db_name}")
    except Exception as e:
        QMessageBox.critical(widget, "Error", f"Failed to load database: {e}")

#   PROT DATABASE
def choose_database(widget, label_to_update):
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QFileDialog, QMessageBox, QProgressDialog, QApplication
    import subprocess
    import os

    file_paths, _ = QFileDialog.getOpenFileNames(
        widget,
        "Choose FASTA Files",
        "",
        "FASTA Files (*.fasta *.fa *.txt);;All Files (*)"
    )

    if not file_paths:
        return

    # Okienko postępu
    progress = QProgressDialog("Building BLAST database...", "Cancel", 0, 0, widget)
    progress.setWindowTitle("Please wait")
    progress.setWindowModality(Qt.WindowModality.ApplicationModal)
    progress.setMinimumDuration(0)
    progress.show()
    QApplication.processEvents()

    merged_path = "merged_database.fasta"

    try:
        with open(merged_path, "w", encoding="utf-8") as outfile:
            for path in file_paths:
                with open(path, "r", encoding="utf-8") as infile:
                    content = infile.read()
                    if not content.endswith("\n"):
                        content += "\n"
                    outfile.write(content)
    except Exception as e:
        progress.close()
        QMessageBox.critical(widget, "Error", f"Failed to merge FASTA files: {e}")
        return

    db_name = os.path.splitext(os.path.basename(merged_path))[0]

    try:
        subprocess.run([
            "BLAST/makeblastdb",
            "-in", merged_path,
            "-dbtype", "prot",
            "-out", f"database/{db_name}"
        ], check=True)

        filenames = [os.path.basename(p) for p in file_paths]
        label_to_update.setText("+\n ".join(filenames))
        QMessageBox.information(widget, "Success", f"Database created.")
        widget.open_blast_button.setEnabled(True)
        widget.blast_button.setEnabled(True)
        widget.blastn_button.setEnabled(True)
        widget.blastx_button.setEnabled(True)
        widget.tblastn_button.setEnabled(True)
        widget.tblastx_button.setEnabled(True)
        
    except Exception as e:
        QMessageBox.critical(widget, "Error", f"Failed to run makeblastdb: {e}")

    finally:
        progress.close()


#   NEUTRIN DATABASE
def choose_database_n(widget, label_to_update):
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QFileDialog, QMessageBox, QProgressDialog, QApplication
    import subprocess
    import os

    file_paths, _ = QFileDialog.getOpenFileNames(
        widget,
        "Choose FASTA Files",
        "",
        "FASTA Files (*.fasta *.fa *.txt);;All Files (*)"
    )

    if not file_paths:
        return

    # Okienko postępu
    progress = QProgressDialog("Building BLAST database...", "Cancel", 0, 0, widget)
    progress.setWindowTitle("Please wait")
    progress.setWindowModality(Qt.WindowModality.ApplicationModal)
    progress.setMinimumDuration(0)
    progress.show()
    QApplication.processEvents()

    merged_path = "merged_database.fasta"

    try:
        with open(merged_path, "w", encoding="utf-8") as outfile:
            for path in file_paths:
                with open(path, "r", encoding="utf-8") as infile:
                    content = infile.read()
                    if not content.endswith("\n"):
                        content += "\n"
                    outfile.write(content)
    except Exception as e:
        progress.close()
        QMessageBox.critical(widget, "Error", f"Failed to merge FASTA files: {e}")
        return

    db_name = os.path.splitext(os.path.basename(merged_path))[0]

    try:
        subprocess.run([
            "BLAST/makeblastdb",
            "-in", merged_path,
            "-dbtype", "nucl",
            "-out", f"database/{db_name}"
        ], check=True)

        filenames = [os.path.basename(p) for p in file_paths]
        label_to_update.setText("+ ".join(filenames))
        QMessageBox.information(widget, "Success", f"Database created.")
        widget.label_database_name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        widget.open_blast_button.setEnabled(True)
        widget.blast_button.setEnabled(True)
        widget.blastn_button.setEnabled(True)
        widget.blastx_button.setEnabled(True)
        widget.tblastn_button.setEnabled(True)
        widget.tblastx_button.setEnabled(True)

    except Exception as e:
        QMessageBox.critical(widget, "Error", f"Failed to run makeblastdb: {e}")

    finally:
        progress.close()
        
        
#/////////////////////////////////////////////////////////////////////////////////////////
########################## RUN BLAST P ###################################################
import tempfile
import os

def run_blast(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox, QProgressDialog, QApplication
    from PyQt6.QtCore import Qt
    import subprocess
    import os
    from blast_parser import parse_blast_output
    from blast_view import parse_blast_xml

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

        #db_name = db_label.text().strip()
       # if db_name.lower() == "none":
          #  QMessageBox.warning(widget, "Error", "No database selected.")
         #   return

        
        progress = QProgressDialog("Running BLAST, please wait...", None, 0, 0, widget)
        progress.setWindowTitle("Running BLAST")
        progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        progress.show()
        QApplication.processEvents()

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/blastp",
            "-query", "s.txt",
            "-db", "database/merged_database",
            "-out", "output.asn1",
            "-outfmt", "11"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.txt",
            "-outfmt", "0"
        ], check=True)

        # PARSE AND SAVE
        parse_blast_output("output.xml", "blast_out.txt")

        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())

        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)

        QMessageBox.information(widget, "BLAST", "BLAST complete!")
        widget.save_sequences_button.setEnabled(True)
        widget.open_blast_file_button.setEnabled(True)

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e))

    finally:
        progress.close()
        
#/////////////////////////////////////////////////////////////////////////////////////////        
########################## RUN BLAST N ###################################################        
def run_blastn(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox, QProgressDialog, QApplication
    from PyQt6.QtCore import Qt
    import subprocess
    import os
    from blast_parser import parse_blast_output
    from blast_view import parse_blast_xml

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

        #db_name = db_label.text().strip()
        #if db_name.lower() == "none":
            #QMessageBox.warning(widget, "Error", "No database selected.")
            #return

        
        progress = QProgressDialog("Running BLAST, please wait...", None, 0, 0, widget)
        progress.setWindowTitle("Running BLAST")
        progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        progress.show()
        QApplication.processEvents()

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/blastn",
            "-query", "s.txt",
            "-db", "database/merged_database",
            "-out", "output.asn1",
            "-outfmt", "11"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.txt",
            "-outfmt", "0"
        ], check=True)

        # PARSE AND SAVE
        parse_blast_output("output.xml", "blast_out.txt")

        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())

        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)

        QMessageBox.information(widget, "BLAST", "BLAST complete!")
        widget.save_sequences_button.setEnabled(True)
        widget.open_blast_file_button.setEnabled(True)

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e))

    finally:
        progress.close()
        
#/////////////////////////////////////////////////////////////////////////////////////////
########################## RUN BLAST X ###################################################
def run_blastx(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox, QProgressDialog, QApplication
    from PyQt6.QtCore import Qt
    import subprocess
    import os
    from blast_parser import parse_blast_output
    from blast_view import parse_blast_xml

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

        #db_name = db_label.text().strip()
        #if db_name.lower() == "none":
        #    QMessageBox.warning(widget, "Error", "No database selected.")
        #    return

        
        progress = QProgressDialog("Running BLAST, please wait...", None, 0, 0, widget)
        progress.setWindowTitle("Running BLAST")
        progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        progress.show()
        QApplication.processEvents()

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/blastx",
            "-query", "s.txt",
            "-db", "database/merged_database",
            "-out", "output.asn1",
            "-outfmt", "11"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.txt",
            "-outfmt", "0"
        ], check=True)

        # PARSE AND SAVE
        parse_blast_output("output.xml", "blast_out.txt")

        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())

        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)

        QMessageBox.information(widget, "BLAST", "BLAST complete!")
        widget.save_sequences_button.setEnabled(True)
        widget.open_blast_file_button.setEnabled(True)

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e))

    finally:
        progress.close()
        
        

def run_tblastn(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox, QProgressDialog, QApplication
    from PyQt6.QtCore import Qt
    import subprocess
    import os
    from blast_parser import parse_blast_output
    from blast_view import parse_blast_xml

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

       # db_name = db_label.text().strip()
       # if db_name.lower() == "none":
        #    QMessageBox.warning(widget, "Error", "No database selected.")
       #     return

        
        progress = QProgressDialog("Running BLAST, please wait...", None, 0, 0, widget)
        progress.setWindowTitle("Running BLAST")
        progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        progress.show()
        QApplication.processEvents()

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/tblastn",
            "-query", "s.txt",
            "-db", "database/merged_database",
            "-out", "output.asn1",
            "-outfmt", "11"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.txt",
            "-outfmt", "0"
        ], check=True)

        # PARSE AND SAVE
        parse_blast_output("output.xml", "blast_out.txt")

        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())

        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)

        QMessageBox.information(widget, "BLAST", "BLAST complete!")
        widget.save_sequences_button.setEnabled(True)
        widget.open_blast_file_button.setEnabled(True)

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e))

    finally:
        progress.close()



def run_tblastx(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox, QProgressDialog, QApplication
    from PyQt6.QtCore import Qt
    import subprocess
    import os
    from blast_parser import parse_blast_output
    from blast_view import parse_blast_xml

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

       # db_name = db_label.text().strip()
       # if db_name.lower() == "none":
       #     QMessageBox.warning(widget, "Error", "No database selected.")
       #     return

        
        progress = QProgressDialog("Running BLAST, please wait...", None, 0, 0, widget)
        progress.setWindowTitle("Running BLAST")
        progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        progress.show()
        QApplication.processEvents()

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/tblastx",
            "-query", "s.txt",
            "-db", "database/merged_database",
            "-out", "output.asn1",
            "-outfmt", "11"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)
        
        subprocess.run([
            "BLAST/blast_formatter",
            "-archive", "output.asn1",
            "-out", "output.txt",
            "-outfmt", "0"
        ], check=True)

        # PARSE AND SAVE
        parse_blast_output("output.xml", "blast_out.txt")

        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())

        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)

        QMessageBox.information(widget, "BLAST", "BLAST complete!")
        widget.save_sequences_button.setEnabled(True)
        widget.open_blast_file_button.setEnabled(True)

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e))

    finally:
        progress.close()

#/////////////////////////////////////////////////////////////////////////////////////////
########################## SAVE SEQUENCES WITH BLAST INFO ################################

def save_sequences_from_blast(widget):
    try:
        from PyQt6.QtWidgets import QFileDialog, QMessageBox

        # Krok 1: Wczytaj sekwencje FASTA jako dict {ID: sequence}
        seq_dict = {}
        with open("merged_database.fasta", "r", encoding="utf-8") as f:
            content = f.read().split(">")
            for entry in content:
                if not entry.strip():
                    continue
                lines = entry.strip().split("\n")
                header = lines[0]
                sequence = "".join(lines[1:])
                seq_id = header.split()[0].strip()  # np. AT1G01010.1
                seq_dict[seq_id] = sequence.strip()

        # Krok 2: Wczytaj ID + pełny wiersz z BLAST
        selected_entries = []  # lista par (id, linia z blast)
        with open("blast_out.txt", "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("Query:") or line.startswith("BLAST Results:") or not line.strip():
                    continue
                identifier = line.split()[0].strip()
                selected_entries.append((identifier, line.strip()))

        # Krok 3: Zapisz sekwencje z nagłówkiem z BLAST
        output_path, _ = QFileDialog.getSaveFileName(widget, "Save FASTA File", "", "FASTA Files (*.fasta *.fa);;All Files (*)")
        if not output_path:
            return

        with open(output_path, "w", encoding="utf-8") as f:
            for seq_id, blast_line in selected_entries:
                sequence = seq_dict.get(seq_id)
                if sequence:
                    f.write(f">{blast_line}\n")
                    f.write(f"{sequence}\n")

        QMessageBox.information(widget, "Saved", "Sequences saved to FASTA file.")

    except Exception as e:
        QMessageBox.critical(widget, "Error", str(e))


import xml.etree.ElementTree as ET


def open_sequence_from_xml(widget,output_textbox):
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
    from blast_parser import parse_blast_output
    import subprocess
    import os

    file_paths, _ = QFileDialog.getOpenFileNames(
        widget,
        "Choose BLAST XML Files",
        "",
        "XML Files (*.xml);;All Files (*)"
    )

    if file_paths:
        first_path = file_paths[0]
        

    if not file_paths:
        return
    
    parse_blast_output(first_path, "blast_out.txt")
    with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())
            
    from blast_view import parse_blast_xml
    query_alignments = parse_blast_xml("output.xml")
    widget.query_list_widget.clear()
    for query in query_alignments:
        widget.query_list_widget.addItem(query)
        
        
def show_blast_file(path):
    import os
    import platform
    import subprocess
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return

    system = platform.system()

    try:
        if system == "Windows":
            os.startfile(path)
        elif system == "Darwin":  # macOS
            subprocess.run(["open", path], check=False)
        elif system == "Linux":
            subprocess.run(["xdg-open", path], check=False)
        else:
            print(f"Unsupported OS: {system}")
    except Exception as e:
        print(f"Failed to open file: {e}")
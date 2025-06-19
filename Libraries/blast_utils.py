from blast_parser import parse_blast_output

import subprocess
from PyQt6.QtWidgets import QFileDialog, QMessageBox


###################### CREATING DATABASE ########################################
#   PROT DATABASE
def choose_database(widget, label_to_update):
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
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
        QMessageBox.critical(widget, "Error", f"Failed to merge FASTA files: {e}")
        return

    db_name = os.path.splitext(os.path.basename(merged_path))[0]

    try:
        subprocess.run([
            "BLAST/makeblastdb",
            "-in", merged_path,
            "-dbtype", "prot",
            "-out", f"baza/{db_name}"
        ], check=True)

        import os

        filenames = [os.path.basename(p) for p in file_paths]
        label_to_update.setText("+ ".join(filenames))
        QMessageBox.information(widget, "Success", f"Database created: {db_name}")

    except Exception as e:
        QMessageBox.critical(widget, "Error", f"Failed to run makeblastdb: {e}")


#   NEUTRIN DATABASE
def choose_database_n(widget, label_to_update):
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
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
        QMessageBox.critical(widget, "Error", f"Failed to merge FASTA files: {e}")
        return

    db_name = os.path.splitext(os.path.basename(merged_path))[0]

    try:
        subprocess.run([
            "BLAST/makeblastdb",
            "-in", merged_path,
            "-dbtype", "nucl",
            "-out", f"baza/{db_name}"
        ], check=True)

        import os

        filenames = [os.path.basename(p) for p in file_paths]
        label_to_update.setText("+ ".join(filenames))
        QMessageBox.information(widget, "Success", f"Database created: {db_name}")

    except Exception as e:
        QMessageBox.critical(widget, "Error", f"Failed to run makeblastdb: {e}")
        
        
#/////////////////////////////////////////////////////////////////////////////////////////
########################## RUN BLAST P ###################################################
import tempfile
import os

def run_blast(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox
    import subprocess
    import os
    from blast_parser import parse_blast_output

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

        db_name = db_label.text().strip()
        if db_name.lower() == "none":
            QMessageBox.warning(widget, "Error", "No database selected.")
            return

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/blastp",
            "-query", "s.txt",
            "-db", f"baza/merged_database",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)

        # READ QUERY FOR MAKING OUTPUT FILE
        seq_dict = {}
        with open("s.txt", "r", encoding="utf-8") as f:
            content = f.read().strip().split(">")
            for entry in content:
                if not entry:
                    continue
                parts = entry.split("\n", 1)
                if len(parts) == 2:
                    header = parts[0].strip()
                    sequence = parts[1].replace("\n", "")
                    seq_dict[header] = sequence

        # PARSE AND SAVE sequences.txt
        parse_blast_output("output.xml", seq_dict, "blast_out.txt")

        # SHOW FILE ON APPLICATION WINDOW
        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())
        
        
        from blast_view import parse_blast_xml
        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)
            
        QMessageBox.information(widget, "BLAST", "BLAST complete!")

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e))
        
#/////////////////////////////////////////////////////////////////////////////////////////        
########################## RUN BLAST N ###################################################        
def run_blastn(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox
    import subprocess
    import os
    from blast_parser import parse_blast_output

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

        db_name = db_label.text().strip()
        if db_name.lower() == "none":
            QMessageBox.warning(widget, "Error", "No database selected.")
            return

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/blastn",
            "-query", "s.txt",
            "-db", f"baza/merged_database",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)

        # READ QUERY FOR MAKING OUTPUT FILE
        seq_dict = {}
        with open("s.txt", "r", encoding="utf-8") as f:
            content = f.read().strip().split(">")
            for entry in content:
                if not entry:
                    continue
                parts = entry.split("\n", 1)
                if len(parts) == 2:
                    header = parts[0].strip()
                    sequence = parts[1].replace("\n", "")
                    seq_dict[header] = sequence

        # PARSE AND SAVE sequences.txt
        parse_blast_output("output.xml", seq_dict, "blast_out.txt")

        # SHOW FILE ON APPLICATION WINDOW
        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())
        
        
        from blast_view import parse_blast_xml
        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)
            
        QMessageBox.information(widget, "BLAST", "BLAST complete!")

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e)) 
        
#/////////////////////////////////////////////////////////////////////////////////////////
########################## RUN BLAST X ###################################################
def run_blastx(widget, input_textbox, output_textbox, db_label):
    from PyQt6.QtWidgets import QMessageBox
    import subprocess
    import os
    from blast_parser import parse_blast_output

    try:
        sequences = input_textbox.toPlainText().strip()
        if not sequences:
            QMessageBox.warning(widget, "Error", "No sequences provided.")
            return

        db_name = db_label.text().strip()
        if db_name.lower() == "none":
            QMessageBox.warning(widget, "Error", "No database selected.")
            return

        # SAVE QUERY FOR BLAST
        with open("s.txt", "w", encoding="utf-8") as f:
            f.write(sequences.replace("\n", os.linesep))

        # RUN BLAST
        subprocess.run([
            "BLAST/blastx",
            "-query", "s.txt",
            "-db", f"baza/merged_database",
            "-out", "output.xml",
            "-outfmt", "5"
        ], check=True)

        # READ QUERY FOR MAKING OUTPUT FILE
        seq_dict = {}
        with open("s.txt", "r", encoding="utf-8") as f:
            content = f.read().strip().split(">")
            for entry in content:
                if not entry:
                    continue
                parts = entry.split("\n", 1)
                if len(parts) == 2:
                    header = parts[0].strip()
                    sequence = parts[1].replace("\n", "")
                    seq_dict[header] = sequence

        # PARSE AND SAVE sequences.txt
        parse_blast_output("output.xml", seq_dict, "blast_out.txt")

        # SHOW FILE ON APPLICATION WINDOW
        with open("blast_out.txt", "r", encoding="utf-8") as f:
            output_textbox.setPlainText(f.read())
        
        
        from blast_view import parse_blast_xml
        query_alignments = parse_blast_xml("output.xml")
        widget.query_list_widget.clear()
        for query in query_alignments:
            widget.query_list_widget.addItem(query)
            
        QMessageBox.information(widget, "BLAST", "BLAST complete!")

    except Exception as e:
        QMessageBox.critical(widget, "BLAST Error", str(e))


#/////////////////////////////////////////////////////////////////////////////////////////
########################## SAVE SEQUENCES WITH BLAST INFO ################################

def save_sequences_from_blast(widget):
    try:
        
        seq_dict = {}
        with open("merged_database.fasta", "r", encoding="utf-8") as f:
            content = f.read().split(">")
            for entry in content:
                if not entry.strip():
                    continue
                lines = entry.strip().split("\n")
                header = lines[0]
                sequence = "".join(lines[1:])
                seq_dict[header.strip()] = sequence.strip()

      
        selected_ids = []
        with open("blast_out.txt", "r", encoding="utf-8") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("Query:") or line.startswith("BLAST Results:") or not line.strip():
                    continue
                identifier = line.split()[0].strip()
                selected_ids.append(line.strip())

       
        output_path, _ = QFileDialog.getSaveFileName(widget, "Save FASTA File", "", "FASTA Files (*.fasta *.fa);;All Files (*)")
        if not output_path:
            return

        with open(output_path, "w", encoding="utf-8") as f:
            for line in selected_ids:
                name = line.split()[0]
                header = ">" + line
                sequence = seq_dict.get(name)
                if sequence:
                    f.write(header + "\n")
                    f.write(sequence + "\n")

        QMessageBox.information(widget, "Saved", "Sequences saved to FASTA file.")

    except Exception as e:
        QMessageBox.critical(widget, "Error", str(e))

import matplotlib.patches
import sys
from PyQt6.QtWidgets import (
    QListWidgetItem,
    QApplication, QMainWindow, QWidget, QTabWidget, QLabel, QVBoxLayout, QHBoxLayout,
    QPushButton, QCheckBox, QListWidget, QTextEdit, QLineEdit, QFileDialog, QSpinBox,
    QStatusBar, QFrame, QSpacerItem, QSizePolicy, QDialog, QVBoxLayout, QLabel, QPlainTextEdit, QPushButton
)
from PyQt6.QtWidgets import QProgressDialog
from PyQt6.QtGui import QPixmap, QIcon, QBrush, QColor
from PyQt6.QtCore import Qt, QItemSelectionModel
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class ChartWidget(QWidget):
    ####################### ------- #####################
    def __init__(self, title):
        super().__init__()
        layout = QVBoxLayout()
        self.label = QLabel(title)
        layout.addWidget(self.label)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.hist_lengths = []
        min_seq_len = 0
        max_seq_len = 0
        min_score = 0
        max_score = 0
        min_eval = 0
        max_eval = 0

############## LIBRARIES #######################
from thread_utils import DeleteDuplicatesWorker
from fasta_utils import load_fasta_file, clean_fasta_entries, save_fasta_to_db, delete_duplicates, fetch_all_sequences, fetch_filtered_sequences, fetch_advanced_filtered_sequences
from draw_utils import draw_length_histogram, draw_bitscore_histogram, draw_evalue_histogram
import matplotlib.patches
from move_utils import move_checked_items
from blast_utils import choose_database, choose_database_n, run_blast,run_blastn,run_blastx, save_sequences_from_blast
from blast_view import parse_blast_xml, BlastAlignmentViewer
 

class MainWindow(QMainWindow):

    ################# STARTUP ###################################
    def __init__(self):
        super().__init__()
        min_seq_len = 0
        max_seq_len = 0
        min_score = 0
        max_score = 0
        min_eval = 0
        max_eval = 0
        self.setWindowTitle("TidyTOOL")
        self.setWindowIcon(QIcon("Resources/tidytool.png"))
        #self.setMinimumSize(1600, 900)

        self.status = QStatusBar()
        self.setStatusBar(self.status)

        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        self.cleaner_tab = QWidget()
        self.blast_tab = QWidget()
        self.tabs.addTab(self.cleaner_tab, "CLEANER")
        self.tabs.addTab(self.blast_tab, "BLAST")

        self.setup_cleaner_tab()
        self.setup_blast_tab()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< UI >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ########## VERTICAL SEPARATOR #################
    def make_separator(self):
        frame = QFrame()
        frame.setFrameShape(QFrame.Shape.VLine)
        frame.setFrameShadow(QFrame.Shadow.Sunken)
        frame.setFixedWidth(2)
        return frame
    ########## HORIZONTAL SEPARATOR #################
    def make_hor_separator(self):
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFrameShadow(QFrame.Shadow.Sunken)
        separator.setFixedHeight(4)
        separator.setStyleSheet("color: gray; background-color: gray;")
        return separator
    ########## HORIZONTAL SPACE #################
    def make_hor_space(self):
        hor_space = QFrame()
        hor_space.setFrameShape(QFrame.Shape.HLine)
        hor_space.setFrameShadow(QFrame.Shadow.Plain)
        hor_space.setFixedHeight(15)
        return hor_space
        
    
    ################
    def update_progress_text(self, text):
        self.progress_dialog.setLabelText(text)
    ################
    def update_progress_text(self, text):
        self.progress_dialog.setLabelText(text)
    ################
    def update_progress_value(self, value):
        self.progress_dialog.setValue(value)
        self.progress_dialog.setLabelText(f"Progress: {value}%")
        
    ##################
    def setup_cleaner_tab(self):
        main_layout = QHBoxLayout()

        ############## Left panel ################
        outer_left_layout = QVBoxLayout()
        outer_left_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        inner_left_widget = QWidget()
        left_layout = QVBoxLayout(inner_left_widget)
        #------------------
        self.label_analyse = QLabel("Analyse FASTA")
        font_bold = self.label_analyse.font()
        font_bold.setBold(True)
        font_bold.setPointSize(14)
        self.label_analyse.setFont(font_bold)
        #------------------
        self.label_file = QLabel("File:")
        self.label_file_name = QLabel("None")
        self.load_button = QPushButton("Load Files")
        self.load_button.setFixedSize(150, 40)
        self.load_button.clicked.connect(self.load_file)
        #------------------
        self.show_files_button = QPushButton("Show Files")
        self.show_files_button.setFixedSize(150, 30)
        self.show_files_button.clicked.connect(self.show_loaded_files)
        self.show_files_button.setEnabled(False)
        #------------------
        self.label_checkdup = QLabel("Check duplicates in")
        font_bold2 = self.label_analyse.font()
        font_bold2.setBold(True)
        font_bold2.setPointSize(12)
        self.label_checkdup.setFont(font_bold2)
        #------------------
        checkbox_row = QHBoxLayout()
        self.name_box = QCheckBox("Names")
        self.name_box.setChecked(True)
        #------------------
        self.sequence_box = QCheckBox("Sequences")
        self.sequence_box.setChecked(True)
        #------------------
        #self.show_removed_checkbox = QCheckBox("Show Removed")
        self.show_removed_button = QPushButton("Show Removed")
        
        self.show_duplicates_button = QPushButton("Show Duplicates")
        self.show_duplicates_button.setFixedSize(150, 30)
        self.show_duplicates_button.clicked.connect(self.show_duplicate_sequences)
        self.show_duplicates_button.setEnabled(False)                                                            
        self.show_removed_button.setFixedSize(150, 30)
        self.show_removed_button.clicked.connect(self.show_removed_sequences)
        self.show_removed_button.setEnabled(False)
        
        checkbox_row.addWidget(self.name_box)
        checkbox_row.addWidget(self.sequence_box)
        checkbox_row.addWidget(self.show_removed_button)
        #------------------
        self.clean_button = QPushButton("Clean")
        self.clean_button.setFixedSize(150, 40)
        self.clean_button.clicked.connect(self.clean_file)
        self.clean_button.setEnabled(False)
        #------------------
        self.label_filter = QLabel("Analyze Parameters")
        self.label_filter.setFont(font_bold2)
        #------------------
        filter_checkbox_row = QHBoxLayout()
        self.protein_checkbox_obj = QCheckBox("Protein")
        self.protein_checkbox_obj.setChecked(True)
        self.protein_checkbox_obj.setEnabled(False)
        self.protein_checkbox_obj.stateChanged.connect(self.update_atg_checkbox_text)
        #------------------
        self.gene_checkbox_obj = QCheckBox("Gene")
        self.gene_checkbox_obj.setEnabled(False)
        self.gene_checkbox_obj.stateChanged.connect(self.update_atg_checkbox_text)
        #------------------
        self.atg_checkbox = QCheckBox("Check M")
        self.atg_checkbox.setEnabled(False)
        #------------------
        filter_checkbox_row.addWidget(self.protein_checkbox_obj)
        filter_checkbox_row.addWidget(self.gene_checkbox_obj)
        filter_checkbox_row.addWidget(self.atg_checkbox)
        #------------------
        self.length_checkbox_obj = QCheckBox("Length")
        self.length_checkbox_obj.setEnabled(False)
        #------------------
        length_row = QHBoxLayout()
        self.len_box_low = QSpinBox()
        self.len_box_low.setMaximum(900000)
        self.len_box_low.setFixedWidth(70)
        self.len_box_hi = QSpinBox()
        self.len_box_hi.setMaximum(900000)
        self.len_box_hi.setFixedWidth(70)
        length_row.addWidget(self.length_checkbox_obj)
        length_row.addWidget(QLabel("MIN:"))
        length_row.addWidget(self.len_box_low)
        length_row.addWidget(QLabel("MAX:"))
        length_row.addWidget(self.len_box_hi)
        #------------------
        
        
        self.score_checkbox_obj = QCheckBox("Score")
        self.score_checkbox_obj.setEnabled(False)
        #------------------
        score_row = QHBoxLayout()
        self.score_box_low = QSpinBox()
        self.score_box_low.setMaximum(900000)
        self.score_box_low.setFixedWidth(70)
        self.score_box_hi = QSpinBox()
        self.score_box_hi.setMaximum(900000)
        self.score_box_hi.setFixedWidth(70)
        score_row.addWidget(self.score_checkbox_obj)
        score_row.addWidget(QLabel("MIN:"))
        score_row.addWidget(self.score_box_low)
        score_row.addWidget(QLabel("MAX:"))
        score_row.addWidget(self.score_box_hi)
        #------------------
        
        self.eval_checkbox_obj = QCheckBox("E-Value")
        self.eval_checkbox_obj.setEnabled(False)
        #------------------
        eval_row = QHBoxLayout()
        self.eval_box_low = QSpinBox()
        self.eval_box_low.setMaximum(900000)
        self.eval_box_low.setFixedWidth(70)
        self.eval_box_hi = QSpinBox()
        self.eval_box_hi.setMaximum(900000)
        self.eval_box_hi.setFixedWidth(70)
        eval_row.addWidget(self.eval_checkbox_obj)
        eval_row.addWidget(QLabel("MIN:"))
        eval_row.addWidget(self.eval_box_low)
        eval_row.addWidget(QLabel("MAX:"))
        eval_row.addWidget(self.eval_box_hi)
        #------------------
        
        

        self.label_min = QLabel("MIN")
        self.label_max = QLabel("MAX")
        #------------------
        self.label_status = QLabel("Status:")
        #------------------
        self.label_status_text = QLabel("No file selected")
        #------------------
        self.label_removed = QLabel("Removed 0 sequences")
        #------------------
        self.label_inout = QLabel("Input: 0 Output: 0")
        #------------------
        self.amount_label = QLabel("Records: 0")
        #------------------
        self.amount2_label = QLabel("Records: 0")
        #------------------
        self.label_copyright = QLabel("© Liszka, Stołowski")
        #------------------
        self.analyze_button = QPushButton("Filter")
        self.analyze_button.clicked.connect(self.filter_all)
        self.analyze_button.setFixedSize(150, 40)
        self.analyze_button.setEnabled(False)
        #------------------
        self.reset_button = QPushButton("RESET")
        self.reset_button.setFixedSize(150, 40)
        self.reset_button.clicked.connect(self.reset_filters)
        self.reset_button.setEnabled(False)
        #------------------
        self.save_fasta_button = QPushButton("Save FASTA")
        self.save_fasta_button.clicked.connect(self.save_selected_to_fasta)
        self.save_fasta_button.setFixedSize(150, 40)
        self.save_fasta_button.setEnabled(False)
        #------------------
        left_layout.addWidget(self.label_analyse, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.load_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.show_files_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        #left_layout.addWidget(self.label_file, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_file_name, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.make_hor_separator())
        left_layout.addWidget(self.label_checkdup, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addLayout(checkbox_row)
        left_layout.addWidget(self.clean_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        #left_layout.addWidget(self.label_status, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_status_text, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_removed, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_inout, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.make_hor_separator())
        left_layout.addWidget(self.label_filter, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addLayout(filter_checkbox_row)
        left_layout.addLayout(length_row)  
        left_layout.addLayout(score_row)  
        left_layout.addLayout(eval_row)  
        left_layout.addWidget(self.make_hor_space())  
        left_layout.addWidget(self.analyze_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.reset_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.save_fasta_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.make_hor_separator())
        #------------------
        self.image_label = QLabel()
        self.image_label.setPixmap(QPixmap("resources/tidytool.png").scaledToWidth(30, Qt.TransformationMode.SmoothTransformation))
        left_layout.addWidget(self.image_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        left_layout.addWidget(self.label_copyright, alignment=Qt.AlignmentFlag.AlignHCenter)
        #------------------
        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding))
        #------------------
        
        outer_left_layout.addWidget(inner_left_widget)
        #------------------
        main_layout.addLayout(outer_left_layout, 1)
        main_layout.addWidget(self.make_separator())
        main_layout.addLayout(self.setup_middle_panel(), 5)
        main_layout.addWidget(self.make_separator())
        main_layout.addLayout(self.setup_right_panel(), 4)
        #------------------
        self.cleaner_tab.setLayout(main_layout)
        
        
     ########## MIDDLE PANEL #################
    def setup_middle_panel(self):
        layout = QVBoxLayout()
        #------------------
        self.label_genes = QLabel("Genes")
        font_bold = self.label_genes.font()
        font_bold.setBold(True)
        font_bold.setPointSize(14)   
        self.label_genes.setFont(font_bold)
        #------------------
        self.genes_list = QListWidget()
        self.genes_list.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        self.genes_list.setStyleSheet("""
        QListWidget::item:selected {
        background: CornflowerBlue;
        color: white;
        }""")
        self.genes_list.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.genes_list.customContextMenuRequested.connect(self.show_sequence_context_menu)
        #------------------
        self.select_all_checkbox = QCheckBox("Select all")
        self.select_all_checkbox.stateChanged.connect(self.toggle_select_all_genes)
        self.label_selected_names = QLabel("Selected names")
        self.selected_names_listbox = QListWidget()
        self.selected_names_listbox.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        self.selected_names_listbox.setStyleSheet("""
        QListWidget::item:selected {
        background: CornflowerBlue;
        color: white;
        }""")
        #------------------
        self.filter_select_all_checkbox = QCheckBox("Select all")
        self.filter_select_all_checkbox.stateChanged.connect(self.toggle_select_all_selected)
        #------------------
        self.label_similarity = QLabel("Similarity [%]")
        #------------------
        self.accuracy_box = QSpinBox()
        self.accuracy_box.setFixedWidth(70)
        self.accuracy_box.setMaximum(100)
        self.accuracy_box.setValue(100)
        #------------------
        self.filter_textbox = QLineEdit()
       # self.filter_textbox.textChanged.connect(self.apply_filters)
        #------------------
        self.sequence_textbox = QLineEdit()
        #self.sequence_textbox.textChanged.connect(self.apply_filters)
        #self.sequence_textbox.textChanged.connect(self.apply_filters)
        #------------------
        self.add_button = QPushButton("ADD")
        self.add_button.clicked.connect(self.add_selected_genes)
        #------------------
        self.remove_button = QPushButton("Remove")
        self.remove_button.clicked.connect(self.remove_selected_items)
        self.reset_button.setFixedSize(150, 40)
        self.reset_button.clicked.connect(self.reset_filters)
        #------------------
        #------------------
        widgets = [
            self.label_genes, self.genes_list, self.amount_label, self.select_all_checkbox, 
            QLabel("Filter by name:"), self.filter_textbox,
            QLabel("Filter by sequence:"), self.sequence_textbox,self.label_similarity, self.accuracy_box, self.add_button, self.make_hor_separator(),
            self.label_selected_names, self.selected_names_listbox, self.amount2_label,self.filter_select_all_checkbox,
            self.remove_button
        ]
        #------------------
        for widget in widgets:
            layout.addWidget(widget)
        return layout
        #------------------
        #------------------
    ############## RIGHT PANEL ##########################
    def setup_right_panel(self):
        layout = QVBoxLayout()
        self.chart1 = ChartWidget("E-Value Histogram")
        self.chart2 = ChartWidget("Score Histogram")
        self.chart3 = ChartWidget("Lengths Histogram")
        #self.chart1.figure.patch.set_facecolor("darkgray")
        #self.chart2.figure.patch.set_facecolor("darkgray")
        #self.chart3.figure.patch.set_facecolor("darkgray")
        layout.addWidget(self.chart3)
        layout.addWidget(self.chart2)
        layout.addWidget(self.chart1)
        return layout
        
        
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
    ######################### LOAD FASTA FILE ########################################
    def load_file(self):
        file_names, _ = QFileDialog.getOpenFileNames(
            self, "Open FASTA Files", "", "FASTA Files (*.fasta *.fa *.txt);;All Files (*)"
        )
        if file_names:
            #self.label_file_name.setText(", ".join([f.split("/")[-1] for f in file_names]))
            self.loaded_files = file_names
            self.label_status_text.setText("Files selected")
            self.clean_button.setEnabled(True)
            self.save_fasta_button.setEnabled(True)
            self.analyze_button.setEnabled(False)
            self.show_files_button.setEnabled(True)

            import os
            import time
            if os.path.exists("sequences.db"):
                os.remove("sequences.db")

            from fasta_utils import load_fasta_file, save_fasta_to_db

            total_files = len(file_names)
            self.progress_dialog = QProgressDialog("Loading sequences...", "Cancel", 0, total_files, self)
            self.progress_dialog.setWindowTitle("Loading FASTA")
            self.progress_dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
            self.progress_dialog.setMinimumDuration(0)
            self.progress_dialog.setValue(0)
            self.progress_dialog.show()
            QApplication.processEvents()
            
            for i, file_path in enumerate(file_names):
                fasta_lines = load_fasta_file(file_path)
                save_fasta_to_db(fasta_lines, "sequences.db", append=True)

                self.progress_dialog.setValue(i + 1)
                QApplication.processEvents()
                time.sleep(0.2)
                if self.progress_dialog.wasCanceled():
                    break

            self.progress_dialog.close()
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.information(self, "Done", "All FASTA files loaded.")

            self.genes_list.clear()
            self.selected_names_listbox.clear()
            self.min_seq_len = float("inf")
            self.max_seq_len = 0
            
            
    def show_loaded_files(self):
        from PyQt6.QtWidgets import QMessageBox
        if hasattr(self, "loaded_files") and self.loaded_files:
            msg = QMessageBox(self)
            msg.setWindowTitle("Loaded FASTA Files")
            msg.setText("Files:\n\n" + "\n".join(self.loaded_files))
            msg.exec()
        else:
            QMessageBox.information(self, "Info", "No files loaded.")
            
            
    ############### RUN CLEANER IN OTHER THREAD #################
    def clean_file(self):
        from PyQt6.QtWidgets import QProgressDialog
        from PyQt6.QtCore import Qt

        check_name = self.name_box.isChecked()
        check_sequence = self.sequence_box.isChecked()

        self.progress_dialog = QProgressDialog("Starting...", "Cancel", 0, 0, self)
        self.progress_dialog.setWindowTitle("Cleaning Database")
        self.progress_dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress_dialog.setMinimumDuration(0)
        self.progress_dialog.show()

        self.worker = DeleteDuplicatesWorker(check_name, check_sequence)
        self.worker.progress_percent.connect(self.update_progress_value)
        self.worker.progress_text.connect(self.update_progress_text)
        self.worker.finished.connect(self.on_cleaning_finished)
        self.worker.start()

    ############# CLEANING END ############################
    def on_cleaning_finished(self, remaining_sequences, total_before, total_after):
        from fasta_utils import fetch_all_sequences                                           
        self.progress_dialog.close()
        self.genes_list.clear()
        self.hist_lengths = []
        self.min_seq_len = float("inf")
        self.max_seq_len = 0

        self.cleaned_sequences = remaining_sequences
        
        from fasta_utils import load_fasta_file
        original_lines = []
        for f in getattr(self, "loaded_files", []):
            original_lines.extend(load_fasta_file(f))
       # Parsuj wszystkie sekwencje bez usuwania duplikatów
        original_sequences = []
        current_header = None
        current_sequence = ""
        for line in original_lines:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_sequence:
                    original_sequences.append((current_header, current_sequence))
                current_header = line
                current_sequence = ""
            else:
                current_sequence += line
        if current_header and current_sequence:
            original_sequences.append((current_header, current_sequence))
            
        cleaned_set = set(seq.upper().replace(" ", "") for _, seq in self.cleaned_sequences)
        self.removed_sequences = [
            (h, s) for h, s in original_sequences
            if s.upper().replace(" ", "") not in cleaned_set
        ]        
                                                    
        for header, sequence in remaining_sequences:
            self.genes_list.addItem(f"{header} [{len(sequence)}]")
            self.hist_lengths.append(len(sequence))
            self.min_seq_len = min(self.min_seq_len, len(sequence))
            self.max_seq_len = max(self.max_seq_len, len(sequence))

        self.label_removed.setText(f"Removed {total_before - total_after} sequences")
        self.label_inout.setText(f"Input: {total_before} Output: {total_after}")
        self.amount_label.setText(f"Records: {total_after}")
        self.label_min.setText(str(self.min_seq_len))
        self.label_max.setText(str(self.max_seq_len))
        self.len_box_low.setValue(self.min_seq_len)
        self.len_box_hi.setValue(self.max_seq_len)

        self.analyze_button.setEnabled(True)
        self.length_checkbox_obj.setEnabled(True)
        self.protein_checkbox_obj.setEnabled(True)
        self.atg_checkbox.setEnabled(True)
        self.gene_checkbox_obj.setEnabled(True)
        self.len_box_hi.setEnabled(True)
        self.len_box_low.setEnabled(True)
        self.save_fasta_button.setEnabled(True)
        self.clean_button.setEnabled(False)
        self.score_checkbox_obj.setEnabled(True)
        self.eval_checkbox_obj.setEnabled(True)
        self.reset_button.setEnabled(True)
        self.show_removed_button.setEnabled(True)
        self.show_duplicates_button.setEnabled(True)                                            

        draw_length_histogram(self.chart3, self.hist_lengths)
        
        headers = [self.genes_list.item(i).text() for i in range(self.genes_list.count())]
        draw_bitscore_histogram(self.chart2, headers)
        draw_evalue_histogram(self.chart1, headers)
        
    def show_removed_sequences(self):
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QListWidget, QPushButton
        import sqlite3

        conn = sqlite3.connect("duplicates.db")
        cursor = conn.cursor()
        cursor.execute("SELECT header, sequence FROM sequences")
        records = cursor.fetchall()
        conn.close()

        #if not records:
            #return

        dialog = QDialog(self)
        dialog.setWindowTitle("Removed Sequences")
        layout = QVBoxLayout(dialog)

        list_widget = QListWidget()
        for header, seq in records:
            list_widget.addItem(f"{header} [{len(seq)}]")

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)

        layout.addWidget(list_widget)
        layout.addWidget(close_btn)

        dialog.setLayout(layout)
        dialog.resize(700, 500)
        dialog.exec()                                                                                  

    def show_duplicate_sequences(self):
        from PyQt6.QtWidgets import QDialog, QVBoxLayout, QListWidget, QPushButton
        import sqlite3

        conn = sqlite3.connect("duplicates.db")
        cursor = conn.cursor()
        cursor.execute("SELECT header, sequence FROM sequences")
        records = cursor.fetchall()
        conn.close()

        dialog = QDialog(self)
        dialog.setWindowTitle("Duplicate Sequences")
        layout = QVBoxLayout(dialog)

        list_widget = QListWidget()
        for header, seq in records:
            list_widget.addItem(f"{header} [{len(seq)}]")

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)

        layout.addWidget(list_widget)
        layout.addWidget(close_btn)

        dialog.setLayout(layout)
        dialog.resize(700, 500)
        dialog.exec()                                                                            
    ############# ADD GENES TO LIST ############################
    def add_selected_genes(self):
        for i in range(self.genes_list.count()):
            item = self.genes_list.item(i)
            if item.isSelected():
                text = item.text()
                if not any(self.selected_names_listbox.item(j).text() == text for j in range(self.selected_names_listbox.count())):
                    self.selected_names_listbox.addItem(text)
        self.amount2_label.setText(f"Records: {self.selected_names_listbox.count()}")
        
        
    ############# SELECT ALL GENES CHECKBOX ############################
    def toggle_select_all_genes(self, state):
        state == Qt.CheckState.Checked
        for i in range(self.genes_list.count()):
            item = self.genes_list.item(i)
            item.setSelected(state)
    ############# SELECT ALL SELECTED CHECKBOX #########################
    def toggle_select_all_selected(self, state):
        state == Qt.CheckState.Checked
        for i in range(self.selected_names_listbox.count()):
            item = self.selected_names_listbox.item(i)
            item.setSelected(state)
    ############## REMOVE SELECTED FROM LIST ############################
    def remove_selected_items(self):
        for i in reversed(range(self.selected_names_listbox.count())):
            item = self.selected_names_listbox.item(i)
            if item.isSelected():
                self.selected_names_listbox.takeItem(i)
        self.amount2_label.setText(f"Records: {self.selected_names_listbox.count()}")
        
    def show_sequence_context_menu(self, position):
        from PyQt6.QtWidgets import QMenu, QMessageBox

        item = self.genes_list.itemAt(position)
        if item is None:
            return

        menu = QMenu()
        show_action = menu.addAction("Show Sequence")
        action = menu.exec(self.genes_list.viewport().mapToGlobal(position))

        if action == show_action:
           
            header_text = item.text().split(" [")[0]
            
            
            from fasta_utils import fetch_all_sequences
            for header, seq in fetch_all_sequences("sequences.db"):
                if header == header_text:
                    self.show_sequence_dialog(header, seq)
                    break

    from PyQt6.QtWidgets import QDialog, QVBoxLayout, QLabel, QPlainTextEdit, QPushButton

    def show_sequence_dialog(self, header, sequence):
        dialog = QDialog(self)
        dialog.setWindowTitle("Selected Sequence")

        layout = QVBoxLayout()

 
        fasta_text = f"{header}\n{sequence}"
        text_edit = QPlainTextEdit()
        text_edit.setPlainText(fasta_text)
        text_edit.setReadOnly(True)
        layout.addWidget(text_edit)

        button_layout = QHBoxLayout()

        close_button = QPushButton("Close")
        close_button.clicked.connect(dialog.accept)
        button_layout.addWidget(close_button)

        layout.addLayout(button_layout)

        dialog.setLayout(layout)
        dialog.resize(700, 400)
        dialog.exec()
    
    ####################### FILTERING ############################
    def apply_filters(self):
        from PyQt6.QtWidgets import QProgressDialog, QApplication
        from fasta_utils import fetch_filtered_sequences

        name_query = self.filter_textbox.text().strip().lower()
        name_terms = [term.strip() for term in name_query.split(",") if term.strip()]
        seq_query = self.sequence_textbox.text().strip().lower()
        seq_terms = [term.strip() for term in seq_query.split(",") if term.strip()]
        threshold = self.accuracy_box.value()

        self.genes_list.clear()

        self.progress_dialog2 = QProgressDialog("Filtering...", "Cancel", 0, 100, self)
        self.progress_dialog2.setWindowTitle("Filtering Sequences")
        self.progress_dialog2.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress_dialog2.setMinimumDuration(0)
        self.progress_dialog2.setValue(0)
        self.progress_dialog2.show()
        QApplication.processEvents()
        
        def update_progress(value):
            self.progress_dialog2.setValue(value)
            QApplication.processEvents()
            if self.progress_dialog2.wasCanceled():
                raise Exception("Filtering cancelled")

        try:
            filtered = fetch_filtered_sequences(
                name_terms=name_terms,
                seq_terms=seq_terms,
                similarity_threshold=threshold,
                db_path="sequences.db",
                progress_callback=update_progress
            )
        except Exception as e:
            self.progress_dialog2.close()
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Aborted", str(e))
            return

        

        for header, sequence in filtered:
            item = QListWidgetItem(f"{header} [{len(sequence)}]")
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            self.genes_list.addItem(item)
        self.progress_dialog2.close()
        self.amount_label.setText(f"Records: {self.genes_list.count()}")
    

    def filter_all(self):
        from PyQt6.QtWidgets import QProgressDialog, QApplication, QMessageBox
        from fasta_utils import fetch_advanced_filtered_sequences

        self.genes_list.clear()

        name_terms = [t.strip().lower() for t in self.filter_textbox.text().strip().split(",") if t.strip()]
        seq_terms = [t.strip().lower() for t in self.sequence_textbox.text().strip().split(",") if t.strip()]
        threshold = self.accuracy_box.value()

        check_length = self.length_checkbox_obj.isChecked()
        min_len = self.len_box_low.value()
        max_len = self.len_box_hi.value()

        check_atg = self.atg_checkbox.isChecked()
        check_m = self.protein_checkbox_obj.isChecked()
        check_atg_dna = self.gene_checkbox_obj.isChecked()

        check_score = getattr(self, "score_checkbox_obj", None)
        if check_score:
            check_score = check_score.isChecked()
            min_score = self.score_box_low.value()
            max_score = self.score_box_hi.value()
        else:
            check_score = False
            min_score = 0
            max_score = 1000

        check_eval = getattr(self, "eval_checkbox_obj", None)
        if check_eval:
            check_eval = check_eval.isChecked()
            min_eval = self.eval_box_low.value()
            max_eval = self.eval_box_hi.value()
        else:
            check_eval = False
            min_eval = 0
            max_eval = 100

        self.progress_dialog = QProgressDialog("Filtering...", "Cancel", 0, 100, self)
        self.progress_dialog.setWindowTitle("Filtering Sequences")
        self.progress_dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
        self.progress_dialog.setMinimumDuration(0)
        self.progress_dialog.setValue(0)
        self.progress_dialog.show()
        QApplication.processEvents()

        def update_progress(val):
            self.progress_dialog.setValue(val)
            QApplication.processEvents()
            if self.progress_dialog.wasCanceled():
                raise Exception("Filtering cancelled")

        try:
            filtered = fetch_advanced_filtered_sequences(
                name_terms=name_terms,
                seq_terms=seq_terms,
                similarity_threshold=threshold,
                check_length=check_length,
                min_len=min_len,
                max_len=max_len,
                check_atg=check_atg,
                check_m=check_m,
                check_atg_dna=check_atg_dna,
                check_score=check_score,
                min_score=min_score,
                max_score=max_score,
                check_eval=check_eval,
                min_eval=min_eval,
                max_eval=max_eval,
                db_path="sequences.db",
                progress_callback=update_progress
            )
        except Exception as e:
            self.progress_dialog.close()
            QMessageBox.warning(self, "Aborted", str(e))
            return

        self.progress_dialog.close()
        self.genes_list.setUpdatesEnabled(False)
        self.genes_list.setSortingEnabled(False)

        for header, sequence in filtered:
            item = QListWidgetItem(f"{header} [{len(sequence)}]")
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            self.genes_list.addItem(item)
        
        self.genes_list.setSortingEnabled(True)
        self.genes_list.setUpdatesEnabled(True)

        self.amount_label.setText(f"Records: {self.genes_list.count()}")
    
    
        
    ################# SAVE FASTA FILE ############################
    def save_selected_to_fasta(self):
        from thread_utils import DeleteDuplicatesWorker
        from fasta_utils import fetch_all_sequences
        file_name, _ = QFileDialog.getSaveFileName(
            self, "Save Selected FASTA", "", "FASTA Files (*.fasta *.fa);;All Files (*)"
        )
        if file_name:
            selected_headers = {self.selected_names_listbox.item(i).text().split(" [")[0]
                                for i in range(self.selected_names_listbox.count())}
            all_sequences = fetch_all_sequences("sequences.db")
            with open(file_name, "w") as f:
                for header, sequence in all_sequences:
                    if header in selected_headers:
                        f.write(f"{header}\n")
                        for i in range(0, len(sequence), 60):
                            f.write(sequence[i:i+60] + "\n")
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.information(self, "Success", "File saved!")
            
            
    ############## CHECKBOXES TEXTS ############################
    def update_atg_checkbox_text(self):
        sender = self.sender()
        if sender == self.protein_checkbox_obj and self.protein_checkbox_obj.isChecked():
            self.gene_checkbox_obj.blockSignals(True)
            self.gene_checkbox_obj.setChecked(False)
            self.gene_checkbox_obj.blockSignals(False)
            self.atg_checkbox.setText("Check M")
        elif sender == self.gene_checkbox_obj and self.gene_checkbox_obj.isChecked():
            self.protein_checkbox_obj.blockSignals(True)
            self.protein_checkbox_obj.setChecked(False)
            self.protein_checkbox_obj.blockSignals(False)
            self.atg_checkbox.setText("Check ATG")
        elif not self.protein_checkbox_obj.isChecked() and not self.gene_checkbox_obj.isChecked():
            self.atg_checkbox.setText("Check")
            
            
    ############## FILTER BY ATG ############################
    def filter_by_atg_mode(self):
        import re
        from fasta_utils import fetch_all_sequences

        self.genes_list.clear()
        name_query = self.filter_textbox.text().strip().lower()
        name_terms = [term.strip() for term in name_query.split(",") if term.strip()]
        seq_query = self.sequence_textbox.text().strip().lower()
        seq_terms = [term.strip() for term in seq_query.split(",") if term.strip()]
        threshold = self.accuracy_box.value()

        check_length = self.length_checkbox_obj.isChecked()
        min_len = self.len_box_low.value()
        max_len = self.len_box_hi.value()

        check_score = self.score_checkbox_obj.isChecked()
        self.min_score = self.score_box_low.value()
        self.max_score = self.score_box_hi.value()

        check_eval = self.eval_checkbox_obj.isChecked()
        self.min_eval = self.eval_box_low.value()
        self.max_eval = self.eval_box_hi.value()

        sequences = fetch_all_sequences("sequences.db")
        count = 0

        for header, sequence in sequences:
            header_lc = header.lower()
            seq_len = len(sequence)

            # NAME
            if name_terms and not all(term in header_lc for term in name_terms):
                continue

            # SEQUENCES
            match_seq = True
            for term in seq_terms:
                term_len = len(term)
                found = False
                for i in range(len(sequence) - term_len + 1):
                    window = sequence[i:i+term_len].lower()
                    match_count = sum(1 for a, b in zip(term, window) if a == b)
                    similarity = (match_count / term_len) * 100
                    if similarity >= threshold:
                        found = True
                        break
                if not found:
                    match_seq = False
                    break
            if not match_seq:
                continue

            # LENGTH
            if check_length and not (min_len <= seq_len <= max_len):
                continue

            # ATG/M
            if self.atg_checkbox.isChecked():
                if self.protein_checkbox_obj.isChecked() and not sequence.startswith("M"):
                    continue
                elif self.gene_checkbox_obj.isChecked() and not sequence.startswith("ATG"):
                    continue

            # Score
            if check_score:
                score_match = re.search(r"Score: ([\d\.]+)", header)
                if not score_match or not (self.min_score <= float(score_match.group(1)) <= self.max_score):
                    continue

            # E-value
            if check_eval:
                eval_match = re.search(r"E-value: ([\d\.eE-]+)", header)
                if not eval_match or not (self.min_eval <= float(eval_match.group(1)) <= self.max_eval):
                    continue

            # ADD TO LIST
            item = QListWidgetItem(f"{header} [{seq_len}]")
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            self.genes_list.addItem(item)
            count += 1

        self.amount_label.setText(f"Records: {count}")
        
        
    ############## RESET FILTERS AND SHOW ALL ############################
    def reset_filters(self):
        self.filter_textbox.clear()
        self.sequence_textbox.clear()
        self.accuracy_box.setValue(100)
        self.len_box_low.setValue(0)
        self.len_box_hi.setValue(900000)
        self.length_checkbox_obj.setChecked(False)
        self.protein_checkbox_obj.setChecked(False)
        self.gene_checkbox_obj.setChecked(False)
        self.atg_checkbox.setChecked(False)
        self.score_checkbox_obj.setChecked(False)
        self.eval_checkbox_obj.setChecked(False)

        self.genes_list.clear()
        self.len_box_low.setValue(self.min_seq_len)
        self.len_box_hi.setValue(self.max_seq_len)
        sequences = fetch_all_sequences("sequences.db")
        for header, sequence in sequences:
            item = QListWidgetItem(f"{header} [{len(sequence)}]")
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            self.genes_list.addItem(item)

        self.amount_label.setText(f"Records: {self.genes_list.count()}")
        
        
    ############## HISTOGRAM CLICK ############################
    def connect_histogram_click(self):
        self.chart3.canvas.mpl_connect("pick_event", self.on_length_bar_click)
        
    ############## HISTOGRAM CLICK ############################
    def on_length_bar_click(self, event):
        if not hasattr(self, "length_bins") or not isinstance(event.artist, matplotlib.patches.Rectangle):
            return
        index = self.length_bar_patches.index(event.artist)
        bin_min = self.length_bins[index]
        bin_max = self.length_bins[index + 1] - 1
        self.len_box_low.setValue(bin_min)
        self.len_box_hi.setValue(bin_max)
    
    def view_selected_alignment(self):
        selected_items = self.query_list_widget.selectedItems()
        if not selected_items:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Warning", "Please select a query.")
            return
        query_name = selected_items[0].text()
        query_alignments = parse_blast_xml("output.xml")
        if query_name not in query_alignments:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Not Found", f"No alignments found for: {query_name}")
            return

        # <- zapamiętaj referencję
        self.viewer_window = BlastAlignmentViewer(query_name, query_alignments[query_name])
        self.viewer_window.show()
        self.viewer_window.raise_()
        self.viewer_window.activateWindow()
        
    ######################## BLAST UI ###############################
    def setup_blast_tab(self):
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        left_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        right_layout = QVBoxLayout()
        right_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        outer_layout = QVBoxLayout()
        outer_layout.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        
        self.blast_title = QLabel("BLAST")
        font_big = self.blast_title.font()
        font_big.setBold(True)
        font_big.setPointSize(20)
        self.blast_title.setFont(font_big)
        outer_layout.addWidget(self.blast_title, alignment=Qt.AlignmentFlag.AlignHCenter)
        outer_layout.addWidget(self.make_hor_separator())
        
        font_bold2 = self.label_analyse.font()
        font_bold2.setBold(True)
        font_bold2.setPointSize(12)
        #---------------------------LEFT---------------------------
        load_base_layout = QHBoxLayout()
        load_base_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        self.load_base_button = QPushButton("LOAD DATABASE")
        self.load_base_button.setFixedSize(150, 40)
        #self.base_select_button.clicked.connect(lambda: load_database(self, self.label_database_name))
        
        load_base_layout.addWidget(self.load_base_button)
        
        
        database_layout = QHBoxLayout()
        database_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.base_select_button = QPushButton("CREATE PROT DATABASE")
        self.base_select_button.setFixedSize(150, 40)
        self.base_select_button.clicked.connect(lambda: choose_database(self, self.label_database_name))
        self.base2_select_button = QPushButton("CREATE GENE DATABASE")
        self.base2_select_button.setFixedSize(150, 40)
        self.base2_select_button.clicked.connect(lambda: choose_database_n(self, self.label_database_name))
        self.label_database_file = QLabel("Database File:")
        self.label_database_file.setFont(font_bold2)
        self.label_database_name = QLabel("none")

        self.label_sequences = QLabel("Query sequences")
        self.label_sequences.setFont(font_bold2)
        self.rich_text_sequences = QTextEdit()
        self.label_result = QLabel("Result")
        self.label_result.setFont(font_bold2)
        self.rich_text_result = QTextEdit()
        self.rich_text_sequences.setMinimumWidth(800)
        self.rich_text_result.setMinimumWidth(800)
        self.rich_text_sequences.setMinimumHeight(150)
        self.rich_text_result.setMinimumHeight(150)


        runblast_layout = QHBoxLayout()
        runblast_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.blast_button = QPushButton("Run BLASTP")
        self.blast_button.setFixedSize(150, 40)
        self.blast_button.clicked.connect(lambda: run_blast(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.blastn_button = QPushButton("Run BLASTN")
        self.blastn_button.setFixedSize(150, 40)
        self.blastn_button.clicked.connect(lambda: run_blastn(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.blastx_button = QPushButton("Run BLASTX")
        self.blastx_button.setFixedSize(150, 40)
        self.blastx_button.clicked.connect(lambda: run_blastx(self, self.rich_text_sequences, self.rich_text_result, self.label_database_name))
        self.open_blast_button = QPushButton("Open BLAST")
        self.save_sequences_button = QPushButton("Save Sequences")
        self.save_sequences_button.setFixedSize(150, 40)
        self.save_sequences_button.clicked.connect(lambda: save_sequences_from_blast(self))
        
        left_layout.addLayout(load_base_layout)
        
        database_layout.addWidget(self.base_select_button)
        database_layout.addWidget(self.base2_select_button)
        
        runblast_layout.addWidget(self.blast_button)
        runblast_layout.addWidget(self.blastn_button)
        runblast_layout.addWidget(self.blastx_button)
        
        widgets = [
            self.label_database_file,
            self.label_database_name, self.label_sequences,
            self.rich_text_sequences,
            
        ]
        left_layout.addLayout(database_layout)
        
        for widget in widgets:
            left_layout.addWidget(widget, alignment=Qt.AlignmentFlag.AlignHCenter)
        
        left_layout.addLayout(runblast_layout)
        
        widgets = [
            
            self.label_result, self.rich_text_result,
            #self.open_blast_button,

            self.save_sequences_button
        ]
        
        for widget in widgets:
            left_layout.addWidget(widget, alignment=Qt.AlignmentFlag.AlignHCenter)
            
        self.query_list_widget = QListWidget()
        self.query_list_widget.setMinimumHeight(150)
        self.view_alignment_button = QPushButton("VIEW")
        self.view_alignment_button.setFixedSize(150, 40)
        self.view_alignment_button.clicked.connect(self.view_selected_alignment)
        right_layout.addWidget(QLabel("Query list"))
        right_layout.addWidget(self.query_list_widget)
        right_layout.addWidget(self.view_alignment_button)
        
        #---------------------------RIGHT---------------------------
        #main_layout.addWidget(self.make_separator())
        main_layout.addLayout(left_layout, 5)
        main_layout.addWidget(self.make_separator())
        main_layout.addLayout(right_layout, 5)
       # main_layout.addWidget(self.make_separator())
        
        
        # --------------------------- QUERY LIST + VIEW BUTTON ---------------------------
       

        outer_layout.addLayout(main_layout)
        self.blast_tab.setLayout(outer_layout)

    
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    #window.resize(1600, 900)
    #self.setFixedSize(1600, 900)
    window.show()
    #window.showMaximized()
    sys.exit(app.exec())


    

    
    
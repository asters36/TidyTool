# This file is part of TidyTool
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the MIT License

from PyQt6.QtCore import QThread, pyqtSignal
import sqlite3

class DeleteDuplicatesWorker(QThread):
    progress_text = pyqtSignal(str)
    progress_percent = pyqtSignal(int)
    finished = pyqtSignal(list, int, int)

    def __init__(self, check_name, check_sequence, db_path="sequences.db"):
        super().__init__()
        self.check_name = check_name
        self.check_sequence = check_sequence
        self.db_path = db_path

    def run(self):
        import time
        import sqlite3
        import gc

        seen = set()
        buffer_duplicates = []
        buffer_size = 1000

        start_time = time.time()
        self.progress_text.emit("Connecting to database...")

        try:
            conn = sqlite3.connect(self.db_path)
            c = conn.cursor()

            c.execute("SELECT COUNT(*) FROM sequences")
            total = c.fetchone()[0]
            self.progress_text.emit(f"Filtering {total} sequences...")

            rows = c.execute("SELECT header, sequence FROM sequences")

            clean_conn = sqlite3.connect("cleaned.db")
            clean_c = clean_conn.cursor()
            clean_c.execute("DROP TABLE IF EXISTS sequences")
            clean_c.execute("CREATE TABLE sequences (id INTEGER PRIMARY KEY AUTOINCREMENT, header TEXT, sequence TEXT)")

            dup_conn = sqlite3.connect("duplicates.db")
            dup_c = dup_conn.cursor()
            dup_c.execute("DROP TABLE IF EXISTS sequences")
            dup_c.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")

            for i, (header, sequence) in enumerate(rows):
                cleaned_header = header
                if self.check_name:
                    if "[Length" in header:
                        cleaned_header = header.split("[Length")[0].strip()

                if self.check_name and self.check_sequence:
                    key = cleaned_header + sequence
                elif self.check_name:
                    key = cleaned_header
                elif self.check_sequence:
                    key = sequence
                else:
                    key = cleaned_header + sequence

                if key not in seen:
                    seen.add(key)
                    clean_c.execute("INSERT INTO sequences (header, sequence) VALUES (?, ?)", (header, sequence))
                else:
                    buffer_duplicates.append((header, sequence))

                    if len(buffer_duplicates) >= buffer_size:
                        dup_c.executemany(
                            "INSERT INTO sequences (header, sequence) VALUES (?, ?)",
                            buffer_duplicates
                        )
                        buffer_duplicates.clear()

                if i % max(1, total // 100) == 0:
                    self.progress_percent.emit(int((i + 1) / total * 100))

            if buffer_duplicates:
                dup_c.executemany(
                    "INSERT INTO sequences (header, sequence) VALUES (?, ?)",
                    buffer_duplicates
                )

            clean_conn.commit()
            dup_conn.commit()
            clean_conn.close()
            dup_conn.close()
            conn.close()

            gc.collect()

            elapsed = time.time() - start_time
            print(f"Cleaning finished in {elapsed:.2f} seconds")

            self.finished.emit([], total, len(seen))

        except Exception as e:
            print("Error in run():", e)
            self.progress_text.emit(f"Error: {e}")
            self.finished.emit([], 0, 0)
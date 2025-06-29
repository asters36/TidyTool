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
        
        seen = {}
        id_map = {}
        removed = []

        start_time = time.time()
        self.progress_text.emit("Connecting to database...")

        try:
            conn = sqlite3.connect(self.db_path)
            c = conn.cursor()

            #Count records
            c.execute("SELECT COUNT(*) FROM sequences")
            total = c.fetchone()[0]
            self.progress_text.emit(f"Filtering {total} sequences...")

            query = "SELECT id, header, sequence FROM sequences"
            rows = c.execute(query)

            for i, (row_id, header, sequence) in enumerate(rows):
                if self.check_name and self.check_sequence:
                    key = header + sequence
                elif self.check_name:
                    key = header
                elif self.check_sequence:
                    key = sequence
                else:
                    key = header + sequence

                if key not in seen:
                    seen[key] = row_id
                    id_map[row_id] = (row_id, header, sequence)
                else:
                    removed.append((header, sequence))

                if i % max(1, total // 100) == 0:
                    self.progress_percent.emit(int((i + 1) / total * 100))

            conn.close()
            
            # Save duplicates
            self.progress_text.emit("Saving duplicates...")
            with sqlite3.connect("duplicates.db") as dup_conn:
                dup_c = dup_conn.cursor()
                dup_c.execute("DROP TABLE IF EXISTS sequences")
                dup_c.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
                if removed:
                    dup_c.executemany("INSERT INTO sequences (header, sequence) VALUES (?, ?)", removed)
                dup_conn.commit()

            del removed
            gc.collect()

            # Save clean database
            self.progress_text.emit("Saving cleaned database...")
            with sqlite3.connect("cleaned.db") as clean_conn:
                clean_c = clean_conn.cursor()
                clean_c.execute("DROP TABLE IF EXISTS sequences")
                clean_c.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
                clean_c.executemany(
                    "INSERT INTO sequences (header, sequence) VALUES (?, ?)",
                    [(h, s) for _, (_, h, s) in id_map.items()]
                )
                clean_conn.commit()

            gc.collect()

            elapsed = time.time() - start_time
            print(f"Cleaning finished in {elapsed:.2f} seconds")

            self.finished.emit([], total, len(seen))

        except Exception as e:
            print("Błąd w run():", e)
            self.progress_text.emit(f"Błąd: {e}")
            self.finished.emit([], 0, 0)
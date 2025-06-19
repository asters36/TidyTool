from PyQt6.QtCore import QThread, pyqtSignal
from fasta_utils import delete_duplicates, fetch_all_sequences
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
        start_time = time.time()

        self.progress_text.emit("Fetching sequences...")

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute("SELECT id, header, sequence FROM sequences")
        rows = c.fetchall()
        total = len(rows)

        self.progress_text.emit(f"Filtering {total} sequences...")
        seen = {}
        removed = []

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
            else:
                removed.append((header, sequence))

            # Emit progress every 1%
            if i % max(1, total // 100) == 0:
                self.progress_percent.emit(int((i + 1) / total * 100))

        to_keep_ids = set(seen.values())
        to_delete = [(row_id,) for (row_id, _, _) in rows if row_id not in to_keep_ids]

        self.progress_text.emit("Deleting duplicates...")
        BATCH_SIZE = 1000
        delete_total = len(to_delete)

        conn.execute("BEGIN TRANSACTION")
        for i in range(0, delete_total, BATCH_SIZE):
            batch = to_delete[i:i + BATCH_SIZE]
            c.executemany("DELETE FROM sequences WHERE id = ?", batch)
            self.progress_percent.emit(100 if delete_total == 0 else int((i + len(batch)) / delete_total * 100))
        conn.commit()

        self.progress_text.emit("Saving duplicates...")
        with sqlite3.connect("duplicates.db") as dup_conn:
            dup_c = dup_conn.cursor()
            dup_c.execute("DROP TABLE IF EXISTS sequences")
            dup_c.execute("CREATE TABLE sequences (header TEXT, sequence TEXT)")
            if removed:
                dup_c.executemany("INSERT INTO sequences (header, sequence) VALUES (?, ?)", removed)
            dup_conn.commit()

        self.progress_text.emit("Loading cleaned data...")
        c.execute("SELECT header, sequence FROM sequences")
        remaining = c.fetchall()
        conn.close()

        elapsed = time.time() - start_time
        print(f"Finished in {elapsed:.2f} seconds")

        self.finished.emit(remaining, total, len(remaining))
# This file is part of TidyTool
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the MIT License

from PyQt6.QtCore import QThread
from filter_worker import FilterWorker
import sqlite3
###################### FILTERING USING 2 THREADS ##############################
def start_parallel_filtering(parent, filter_args, on_done_callback):
    try:
        conn = sqlite3.connect("cleaned.db")
        c = conn.cursor()
        c.execute("SELECT MIN(rowid), MAX(rowid) FROM sequences")
        min_id, max_id = c.fetchone()
        conn.close()

        if min_id is None or max_id is None:
            on_done_callback([])
            return

        mid_id = (min_id + max_id) // 2
        id_ranges = [(min_id, mid_id), (mid_id + 1, max_id)]

        parent.threads = []
        parent.workers = []
        parent.partial_results = []
        parent.finished_threads = 0

        def handle_result(result):
            parent.partial_results.extend(result)
            parent.finished_threads += 1
            if parent.finished_threads == 2:
                on_done_callback(parent.partial_results)

        for id_start, id_end in id_ranges:
            thread = QThread()
            worker = FilterWorker((id_start, id_end), filter_args)
            worker.moveToThread(thread)

            worker.finished.connect(handle_result)
            thread.started.connect(worker.run)
            worker.finished.connect(thread.quit)
            worker.finished.connect(worker.deleteLater)
            thread.finished.connect(thread.deleteLater)

            parent.threads.append(thread)
            parent.workers.append(worker)
            thread.start()

    except Exception as e:
        print("Error in start_parallel_filtering:", e)
        on_done_callback([])
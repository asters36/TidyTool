# This file is part of TidyTool
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the MIT License

from PyQt6.QtCore import QThread, pyqtSignal
import sqlite3
import re

class GeneLoaderWorker(QThread):
    finished = pyqtSignal(list, list, list, dict)  # headers, lengths, al_lengths, stats

    def run(self):
        headers = []
        lengths = []
        al_lengths = []
        stats = {
            "scores": [],
            "e_values": [],
            "identities": [],
            "positives": [],
            "min_len": float("inf"),
            "max_len": 0
        }

        try:
            with sqlite3.connect("cleaned.db") as conn:
                c = conn.cursor()
                c.execute("SELECT header, sequence FROM sequences")

                for header, sequence in c:
                    seq_len = len(sequence)
                    lengths.append(seq_len)
                    headers.append(f"{header} [{seq_len}]")
                    stats["min_len"] = min(stats["min_len"], seq_len)
                    stats["max_len"] = max(stats["max_len"], seq_len)

                    header_lower = header.lower()

                    def extract(pattern, target, cast_func):
                        m = re.search(pattern, target)
                        if m:
                            try:
                                return cast_func(m.group(1))
                            except:
                                return None
                        return None

                    score = extract(r"score[:=]?\s*([\d.]+)", header_lower, int)
                    if score is not None:
                        stats["scores"].append(score)

                    e_val = extract(r"e[- ]?value[:=]?\s*([\d.eE+-]+)", header_lower, float)
                    if e_val is not None:
                        stats["e_values"].append(e_val)

                    aln_len = extract(r"alignment length[:=]?\s*(\d+)", header_lower, int)
                    if aln_len is not None:
                        al_lengths.append(aln_len)

                    idt = extract(r"identities[:=]?\s*(\d+)", header_lower, int)
                    if idt is not None:
                        stats["identities"].append(idt)

                    pos = extract(r"positives[:=]?\s*(\d+)", header_lower, int)
                    if pos is not None:
                        stats["positives"].append(pos)

        except Exception as e:
            print("GeneLoaderWorker error:", e)

        self.finished.emit(headers, lengths, al_lengths, stats)
# This file is part of TidyTool
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the MIT License

import sqlite3

def load_fasta_file(filepath):
    with open(filepath, "r") as f:
        return f.readlines()

def make_key(header, sequence, check_name, check_sequence):
    if check_name and check_sequence:
        return f"{header}\n{sequence}"
    elif check_name:
        return header
    elif check_sequence:
        return sequence
    return f"{header}\n{sequence}"

###################### CREATING DATABASE ##############################    
def save_fasta_to_db(fasta_path, db_path="sequences.db", append=False):
    with sqlite3.connect(db_path) as conn, open(fasta_path, "r") as f:
        c = conn.cursor()

        if not append:
            c.execute("DROP TABLE IF EXISTS sequences")
            c.execute("""
                CREATE TABLE sequences (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    header TEXT,
                    sequence TEXT
                )
            """)
        else:
            c.execute("""
                CREATE TABLE IF NOT EXISTS sequences (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    header TEXT,
                    sequence TEXT
                )
            """)

        current_header = None
        current_sequence = []

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header and current_sequence:
                    c.execute("INSERT INTO sequences (header, sequence) VALUES (?, ?)", 
                              (current_header[1:], ''.join(current_sequence)))
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

        
        if current_header and current_sequence:
            c.execute("INSERT INTO sequences (header, sequence) VALUES (?, ?)", 
                      (current_header[1:], ''.join(current_sequence)))

        conn.commit()

###################### GET SEQUENCES ############################## //NOT USED
def fetch_all_sequences(db_path="sequences.db"):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("SELECT header, sequence FROM sequences")
    rows = c.fetchall()
    conn.close()
    return rows

###################### FILTERING ##############################
def fetch_advanced_filtered_sequences(                
    name_terms=None,
    seq_terms=None,
    similarity_threshold=100,
    check_length=False,
    min_len=0,
    max_len=1_000_000,
    check_atg=False,
    check_m=False,
    check_atg_dna=False,
    check_score=False,
    min_score=0.0,
    max_score=10000.0,
    check_eval=False,
    min_eval=0.0,
    max_eval=100.0,
    check_alength=False,
    min_alength=0,
    max_alength=1000000,
    check_identities=False,
    min_identities=0,
    max_identities=1000000,
    check_positives=False,
    min_positives=0,
    max_positives=1000000,
    db_path="cleaned.db",
    progress_callback=None
):
    import sqlite3
    import re

    result = []

    parsed_name_terms = []
    for line in name_terms or []:
        terms = [t.strip().lower() for t in line.split(",") if t.strip()]
        if terms:
            parsed_name_terms.append(terms)

    def header_matches(header):
        header_l = header.lower()
        for and_terms in parsed_name_terms:
            if all(term in header_l for term in and_terms):
                return True
        return not parsed_name_terms  # jeśli brak zapytań, wszystko przepuszcza

 
    seq_term_data = [(term.lower(), len(term)) for term in seq_terms or []]

    with sqlite3.connect(db_path) as conn:
        c = conn.cursor()
        c2 = conn.cursor()

        c2.execute("SELECT COUNT(*) FROM sequences")
        total = c2.fetchone()[0]
        processed = 0

        for header, sequence in c.execute("SELECT header, sequence FROM sequences"):
            processed += 1

            if check_length:
                seq_len = len(sequence)
                if seq_len < min_len or seq_len > max_len:
                    continue

            if not header_matches(header):
                continue

            if seq_term_data:
                seq_lower = sequence.lower()
                match_seq = False
                for term, term_len in seq_term_data:
                    if len(seq_lower) < term_len:
                        continue
                    max_mismatches = int((1 - similarity_threshold / 100) * term_len)
                    for i in range(len(seq_lower) - term_len + 1):
                        window = seq_lower[i:i+term_len]
                        mismatches = sum(1 for a, b in zip(term, window) if a != b)
                        if mismatches <= max_mismatches:
                            match_seq = True
                            break
                    if match_seq:
                        break
                if not match_seq:
                    continue

            # ATG / M filtry
            if check_atg:
                if check_m and not sequence.startswith("M"):
                    continue
                if check_atg_dna and not sequence.startswith("ATG"):
                    continue

            # Score
            if check_score:
                match = re.search(r"Score: ([\d\.]+)", header)
                if not match or not (min_score <= float(match.group(1)) <= max_score):
                    continue

            # E-value
            if check_eval:
                match = re.search(r"E-?value: ([\d\.eE+-]+)", header)
                if not match or not (min_eval <= float(match.group(1)) <= max_eval):
                    continue

            # Alignment Length
            if check_alength:
                match = re.search(r"Alignment Length: (\d+)", header)
                if not match or not (min_alength <= int(match.group(1)) <= max_alength):
                    continue

            # Identities
            if check_identities:
                match = re.search(r"Identities: (\d+)", header)
                if not match or not (min_identities <= int(match.group(1)) <= max_identities):
                    continue

            # Positives
            if check_positives:
                match = re.search(r"Positives: (\d+)", header)
                if not match or not (min_positives <= int(match.group(1)) <= max_positives):
                    continue

            result.append((header, sequence))

            if progress_callback and processed % max(1, total // 100) == 0:
                progress_callback(int(processed / total * 100))

        if progress_callback:
            progress_callback(100)

    return result
    

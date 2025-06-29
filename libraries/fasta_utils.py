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




def save_fasta_to_db(fasta_lines, db_path="sequences.db", append=False):
    import sqlite3

    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    if not append:
        c.execute("DROP TABLE IF EXISTS sequences")
        c.execute("CREATE TABLE sequences (id INTEGER PRIMARY KEY AUTOINCREMENT, header TEXT, sequence TEXT)")
    else:
        c.execute("CREATE TABLE IF NOT EXISTS sequences (id INTEGER PRIMARY KEY AUTOINCREMENT, header TEXT, sequence TEXT)")

    current_header = None
    current_sequence = ""

    for line in fasta_lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header and current_sequence:
                c.execute("INSERT INTO sequences (header, sequence) VALUES (?, ?)", (current_header, current_sequence))
            current_header = line
            current_sequence = ""
        else:
            current_sequence += line

    if current_header and current_sequence:
        c.execute("INSERT INTO sequences (header, sequence) VALUES (?, ?)", (current_header, current_sequence))

    conn.commit()
    conn.close()

def fetch_all_sequences(db_path="sequences.db"):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("SELECT header, sequence FROM sequences")
    rows = c.fetchall()
    conn.close()
    return rows


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

    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    
    # SQL-level filtering: header + length
    query = "SELECT header, sequence FROM sequences"
    conditions = []
    params = []

    if name_terms:
        or_blocks = []
        for line in name_terms:
            and_terms = [t.strip().lower() for t in line.split(",") if t.strip()]
            and_conditions = []
            for term in and_terms:
                and_conditions.append("LOWER(header) LIKE ?")
                params.append(f"%{term}%")
            if and_conditions:
                or_blocks.append("(" + " AND ".join(and_conditions) + ")")
        if or_blocks:
            conditions.append("(" + " OR ".join(or_blocks) + ")")


    if check_length:
        conditions.append("LENGTH(sequence) BETWEEN ? AND ?")
        params.extend([min_len, max_len])

    if conditions:
        query += " WHERE " + " AND ".join(conditions)

    c.execute(query, params)
    rows = c.fetchall()

    if not seq_terms and not check_atg and not check_score and not check_eval and not check_alength and not check_identities and not check_positives:
        if progress_callback:
            progress_callback(100)
        return rows

    seq_term_data = [(term.lower(), len(term)) for term in seq_terms or []]
    result = []

    total = len(rows)
    for idx, (header, sequence) in enumerate(rows):
        seq_lower = sequence.lower()
        # SEQUENCE similarity
        if seq_term_data:
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
       
        # ATG / M filters
        if check_atg:
            if check_m and not sequence.startswith("M"):
                continue
            if check_atg_dna and not sequence.startswith("ATG"):
                continue

        # SCORE (parsed from header)
        if check_score:
            score_match = re.search(r"Score: ([\d\.]+)", header)
            if not score_match or not (min_score <= float(score_match.group(1)) <= max_score):
                continue

        # E-VALUE (parsed from header)
        if check_eval:
            eval_match = re.search(r"E-value: ([\d\.eE-]+)", header)
            if not eval_match or not (min_eval <= float(eval_match.group(1)) <= max_eval):
                continue
        
        if check_alength:
            al_match = re.search(r"Alignment Length: (\d+)", header)
            if not al_match or not (min_alength <= int(al_match.group(1)) <= max_alength):
                continue
                
        if check_identities:
            id_match = re.search(r"Identities: (\d+)", header)
            if not id_match or not (min_identities <= int(id_match.group(1)) <= max_identities):
                continue
                
        if check_positives:
            pos_match = re.search(r"Positives: (\d+)", header)
            if not pos_match or not (min_positives <= int(pos_match.group(1)) <= max_positives):
                continue

        result.append((header, sequence))

        if progress_callback and idx % max(1, total // 100) == 0:
            progress_callback(int((idx + 1) / total * 100))

    return result
    

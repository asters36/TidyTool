def load_fasta_file(filepath):
    with open(filepath, "r") as f:
        return f.readlines()

def clean_fasta_entries(lines, check_name, check_sequence):
    entries = []
    unique = set()
    current_header = None
    current_sequence = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header and current_sequence:
                key = make_key(current_header, current_sequence, check_name, check_sequence)
                if key not in unique:
                    unique.add(key)
                    entries.append((current_header, current_sequence))
            current_header = line
            current_sequence = ""
        else:
            current_sequence += line
    if current_header and current_sequence:
        key = make_key(current_header, current_sequence, check_name, check_sequence)
        if key not in unique:
            unique.add(key)
            entries.append((current_header, current_sequence))
    return entries

def make_key(header, sequence, check_name, check_sequence):
    if check_name and check_sequence:
        return f"{header}\n{sequence}"
    elif check_name:
        return header
    elif check_sequence:
        return sequence
    return f"{header}\n{sequence}"


import sqlite3

def save_fasta_to_db(fasta_lines, db_path="sequences.db", append=False):
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

import sqlite3

def fetch_filtered_sequences(name_terms=None, seq_terms=None, similarity_threshold=100,
                             db_path="sequences.db", progress_callback=None):
    import sqlite3
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    query = "SELECT header, sequence FROM sequences"
    conditions = []
    params = []

    if name_terms:
        for term in name_terms:
            conditions.append("LOWER(header) LIKE ?")
            params.append(f"%{term.lower()}%")

    if conditions:
        query += " WHERE " + " AND ".join(conditions)

    c.execute(query, params)
    rows = c.fetchall()

    if not seq_terms:
        return rows

    filtered = []
    prepped_terms = [(term.lower(), len(term)) for term in seq_terms if term]

    total = len(rows)
    for idx, (header, sequence) in enumerate(rows):
        seq_lower = sequence.lower()
        for term, term_len in prepped_terms:
            if len(seq_lower) < term_len:
                continue
            max_mismatches = int((1 - similarity_threshold / 100) * term_len)

            for i in range(len(seq_lower) - term_len + 1):
                window = seq_lower[i:i+term_len]
                mismatches = sum(1 for a, b in zip(term, window) if a != b)
                if mismatches <= max_mismatches:
                    filtered.append((header, sequence))
                    break
            else:
                continue
            break

        if progress_callback and idx % max(1, total // 100) == 0:
            progress_callback(int((idx + 1) / total * 100))

    return filtered
    
    
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
    db_path="sequences.db",
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
        for term in name_terms:
            conditions.append("LOWER(header) LIKE ?")
            params.append(f"%{term.lower()}%")

    if check_length:
        conditions.append("LENGTH(sequence) BETWEEN ? AND ?")
        params.extend([min_len, max_len])

    if conditions:
        query += " WHERE " + " AND ".join(conditions)

    c.execute(query, params)
    rows = c.fetchall()

    if not seq_terms and not check_atg and not check_score and not check_eval:
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

        result.append((header, sequence))

        if progress_callback and idx % max(1, total // 100) == 0:
            progress_callback(int((idx + 1) / total * 100))

    return result
    
    
    
    

def delete_duplicates(check_name=True, check_sequence=True, db_path="sequences.db"):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("SELECT id, header, sequence FROM sequences")
    rows = c.fetchall()

    seen = set()
    to_keep = []
    for row in rows:
        row_id, header, sequence = row
        key = ""
        if check_name and check_sequence:
            key = header + sequence
        elif check_name:
            key = header
        elif check_sequence:
            key = sequence
        else:
            key = header + sequence

        if key not in seen:
            seen.add(key)
            to_keep.append(row_id)

    if to_keep:
        placeholders = ",".join("?" for _ in to_keep)
        c.execute(f"DELETE FROM sequences WHERE id NOT IN ({placeholders})", to_keep)

    conn.commit()
    conn.close()


def delete_duplicates(check_name=True, check_sequence=True, db_path="sequences.db"):
    import sqlite3
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Ustal kolumny unikalności
    if check_name and check_sequence:
        cols = "header, sequence"
    elif check_name:
        cols = "header"
    elif check_sequence:
        cols = "sequence"
    else:
        cols = "header, sequence"

    # Utwórz tymczasową tabelę bez duplikatów
    c.execute("DROP TABLE IF EXISTS sequences_unique")
    c.execute(f"""
        CREATE TABLE sequences_unique AS
        SELECT MIN(id) as id, header, sequence
        FROM sequences
        GROUP BY {cols}
    """)

    # Wyczyść oryginalną tabelę i skopiuj z unikalnej
    c.execute("DELETE FROM sequences")
    c.execute("INSERT INTO sequences (id, header, sequence) SELECT id, header, sequence FROM sequences_unique")
    c.execute("DROP TABLE sequences_unique")

    conn.commit()
    conn.close()


# === OPTIMIZED DUPLICATE REMOVAL ===

def remove_duplicate_sequences(fasta_lines):
    """
    Usuwa duplikaty sekwencji z danych FASTA.
    Zwraca listę unikalnych (header, sequence).
    """
    from itertools import groupby

    records = []
    current_header = None
    current_sequence = []

    for line in fasta_lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                records.append((current_header, "".join(current_sequence)))
            current_header = line
            current_sequence = []
        else:
            current_sequence.append(line)

    if current_header:
        records.append((current_header, "".join(current_sequence)))

    seen = set()
    unique_records = []
    for header, seq in records:
        key = seq.upper().replace(" ", "")  # normalizacja
        if key not in seen:
            seen.add(key)
            unique_records.append((header, seq))

    return unique_records


# === SQLITE-BASED FASTA CLEANING ===

import sqlite3

def clean_fasta_with_sqlite(db_path, min_len=0, max_len=100000):
    """
    Czyści bazę SQLite z duplikatów i filtruje sekwencje po długości.
    Zwraca listę unikalnych (header, sequence).
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Tabela tymczasowa na unikalne sekwencje
    cursor.execute("DROP TABLE IF EXISTS unique_sequences")
    cursor.execute("""
        CREATE TABLE unique_sequences (
            header TEXT,
            sequence TEXT UNIQUE
        )
    """)

    # Wstawiamy tylko unikalne sekwencje spełniające warunki długości
    cursor.execute("SELECT header, sequence FROM sequences")
    count = 0
    for header, sequence in cursor.fetchall():
        seqlen = len(sequence)
        if min_len <= seqlen <= max_len:
            try:
                cursor.execute("INSERT INTO unique_sequences (header, sequence) VALUES (?, ?)", (header, sequence))
                count += 1
            except sqlite3.IntegrityError:
                continue  # Duplikat

    conn.commit()

    # Odczytujemy wynik
    cursor.execute("SELECT header, sequence FROM unique_sequences")
    result = cursor.fetchall()

    conn.close()
    return result

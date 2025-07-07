# This file is part of BLASTnBRUSH
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the GPL v3.0 License

import xml.etree.ElementTree as ET

def parse_blast_output(blast_file_path, output_txt_path):
    tree = ET.parse(blast_file_path)
    root = tree.getroot()

    with open(output_txt_path, "w", encoding="utf-8") as out:
        for iteration in root.findall("BlastOutput_iterations/Iteration"):
            query_def = iteration.findtext("Iteration_query-def")
            out.write(f"Query: {query_def}\n")
            out.write("BLAST Results:\n")

            results = []

            for hit in iteration.findall("Iteration_hits/Hit"):
                hit_def = hit.findtext("Hit_def")
                hit_len = hit.findtext("Hit_len")
                for hsp in hit.findall("Hit_hsps/Hsp"):
                    score = hsp.findtext("Hsp_score", "?")
                    evalue = hsp.findtext("Hsp_evalue", "?")
                    align_len = hsp.findtext("Hsp_align-len", "?")
                    identities = hsp.findtext("Hsp_identity", "?")
                    positives = hsp.findtext("Hsp_positive", "?")
                    gaps = hsp.findtext("Hsp_gaps", "?")

                    try:
                        align_len_val = int(align_len)
                        identities_pct = f"{(int(identities) / align_len_val * 100):.0f}%" if align_len_val else "?"
                        positives_pct = f"{(int(positives) / align_len_val * 100):.0f}%" if align_len_val else "?"
                    except (ValueError, TypeError):
                        identities_pct = "?"
                        positives_pct = "?"

                    results.append((
                        hit_def, hit_len, score, evalue, align_len,
                        identities_pct, positives_pct, gaps
                    ))

            for r in results:
                tag = (
                    f"[Length: {r[1]}, Score: {r[2]}, E-Value: {r[3]}, "
                    f"Alignment Length: {r[4]}, Identities: {r[5]}, Positives: {r[6]}, Gaps: {r[7]}]"
                )
                out.write(f"{r[0]} {tag}\n")

            out.write("\n")
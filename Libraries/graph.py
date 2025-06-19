from Bio.Blast import NCBIXML
import matplotlib.pyplot as plt

# Wczytaj plik BLAST XML
with open("output.xml") as handle:
    blast_records = list(NCBIXML.parse(handle))

# Zbieranie trafień
data = []
for record in blast_records:
    for alignment in record.alignments:
        best_hsp = max(alignment.hsps, key=lambda h: h.score)
        data.append((alignment.hit_def, best_hsp.score, best_hsp.expect, best_hsp.identities))

# Sortowanie po wyniku i wybór top 10
data.sort(key=lambda x: -x[1])
top_data = data[:10]

# Wydzielenie danych
names = [d[0] for d in top_data]
scores = [d[1] for d in top_data]
evalues = [d[2] for d in top_data]
identities = [d[3] for d in top_data]

# Wykresy
plt.figure(figsize=(10, 6))
plt.barh(names, scores)
plt.xlabel("BLAST Score")
plt.title("Top 10 Hits by Score")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.barh(names, evalues)
plt.xlabel("E-value (log scale)")
plt.xscale("log")
plt.title("Top 10 Hits by E-value")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.barh(names, identities)
plt.xlabel("Number of Identical Residues")
plt.title("Top 10 Hits by Identity")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()
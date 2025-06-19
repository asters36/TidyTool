
import xml.etree.ElementTree as ET
from collections import defaultdict
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QScrollArea, QSizePolicy

def parse_blast_xml(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    iterations = root.findall("BlastOutput_iterations/Iteration")

    query_alignments = {}

    for iteration in iterations:
        query_def = iteration.findtext("Iteration_query-def")
        query_len = int(iteration.findtext("Iteration_query-len"))
        alignment_blocks = []

        for hit in iteration.findall("Iteration_hits/Hit"):
            hit_def = hit.findtext("Hit_def")
            for hsp in hit.findall("Hit_hsps/Hsp"):
                query_from = int(hsp.findtext("Hsp_query-from"))
                query_to = int(hsp.findtext("Hsp_query-to"))
                aln_length = abs(query_to - query_from)
                alignment_blocks.append({
                    "query_len": int(query_len),
                    "hit_def": hit_def,
                    "from": query_from,
                    "to": query_to,
                    "alignment_length": aln_length
                })
        query_alignments[query_def] = alignment_blocks

    return query_alignments

class BlastAlignmentViewer(QWidget):
    def __init__(self, query_def, alignments, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Alignments: {query_def}")
        self.showMaximized()

        layout = QVBoxLayout(self)

        grouped = defaultdict(list)
        for aln in alignments:
            grouped[aln["hit_def"]].append(aln)

        # Sortujemy wg maksymalnej długości dopasowania
        sorted_hits = sorted(grouped.items(), key=lambda item: -max(a["alignment_length"] for a in item[1]))

        fig_height = max(4, len(sorted_hits) * 0.4)
        self.canvas = FigureCanvas(Figure(figsize=(12, fig_height)))
        self.toolbar = NavigationToolbar(self.canvas, self)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)

        plot_container = QWidget()
        plot_layout = QVBoxLayout(plot_container)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)

        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.canvas.setMinimumHeight(int(self.canvas.figure.get_size_inches()[1] * self.canvas.figure.dpi))

        scroll_area.setWidget(plot_container)
        layout.addWidget(scroll_area)

        self.setLayout(layout)
        self.plot_alignments(query_def, sorted_hits)

    def plot_alignments(self, query_def, sorted_alignments):
        fig = self.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)

        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(query_def)

        all_alignments = [aln for _, alns in sorted_alignments for aln in alns]
        query_len = max((aln["query_len"] for aln in all_alignments), default=100)
        ax.set_xticks(range(0, query_len + 1, max(10, query_len // 20)))
        ax.set_xticklabels(range(0, query_len + 1, max(10, query_len // 20)), fontsize=8)

        bar_height = 0.12
        row_spacing = 0.5
        text_offset = 0.2
        total = len(sorted_alignments)

        for i, (hit_def, alignments) in enumerate(sorted_alignments):
            y = (total - 1 - i) * row_spacing  # 👈 odwracamy kolejność
            for aln in alignments:
                start = min(aln["from"], aln["to"])
                length = abs(aln["to"] - aln["from"])
                ax.broken_barh([(start, length)], (y - bar_height / 2, bar_height), facecolors='blue')
            ax.text(-5, y, hit_def, ha='right', va='center', fontsize=8, color='black')

        ax.set_xlim(0, query_len)
        ax.set_ylim(0, total * row_spacing + 0.5)
        ax.set_yticks([])
        ax.grid(True, linestyle=':', linewidth=0.5)

        fig.tight_layout(pad=0.5)
        self.canvas.draw()

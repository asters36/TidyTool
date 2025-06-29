# This file is part of TidyTool
# Copyright (c) 2025 Aleksandra Liszka, Artur Stołowski, Aleksandra Marcisz
# Licensed under the MIT License

def draw_length_histogram(chart_widget, lengths, bin_count=5):
    if not lengths or len(lengths) < 2:
        return

    min_len = min(lengths)
    max_len = max(lengths)
    bin_size = max(1, round((max_len - min_len) / bin_count))

    bins = list(range(min_len, max_len + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for length in lengths:
        for i in range(len(bins) - 1):
            if bins[i] <= length < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if length == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bar_container = ax.bar(labels, counts, picker=True)
    ax.set_ylim(0, max(counts) * 1.15)
    for rect, count in zip(bar_container, counts):
        height = rect.get_height()
        if height > 0:
            ax.text(rect.get_x() + rect.get_width() / 2, height + 1,
                    str(count), ha='center', va='bottom', fontsize=8)
    
    ax.set_title("Length")
    ax.set_xlabel("Length range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()

    parent = chart_widget.parentWidget()
    if parent:
        parent.length_bins = bins
        parent.length_bar_patches = list(bar_container.patches)


def draw_bitscore_histogram(chart_widget, headers):
    import re

    scores = []
    for header in headers:
        match = re.search(r"Score: ([\d\.]+)", header)
        if match:
            try:
                scores.append(float(match.group(1)))
            except ValueError:
                pass

    if not scores or len(scores) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Scores Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_score = min(scores)
    max_score = max(scores)
    bin_count = 5
    bin_size = max(1, round((max_score - min_score) / bin_count))
    bins = list(range(int(min_score), int(max_score) + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for score in scores:
        for i in range(len(bins) - 1):
            if bins[i] <= score < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if score == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bars = ax.bar(labels, counts)
    ax.set_ylim(0, max(counts) * 1.15)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    height + 1,
                    str(count),
                    ha='center',
                    va='bottom',
                    fontsize=8)
    
                    
    ax.set_title("Bit Score")
    ax.set_xlabel("Score range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()
    
    
def draw_evalue_histogram(chart_widget, headers):
    import re

    evalues = []
    for header in headers:
        match = re.search(r"E-Value: ([\d\.]+)", header)
        if match:
            try:
                evalues.append(float(match.group(1)))
            except ValueError:
                pass

    if not evalues or len(evalues) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No E-Values Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return
        
    min_score = min(evalues)
    max_score = max(evalues)
    bin_count = 5
    bin_size = max(1, round((max_score - min_score) / bin_count))
    bins = list(range(int(min_score), int(max_score) + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for evalue in evalues:
        for i in range(len(bins) - 1):
            if bins[i] <= evalue < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if evalue == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bars = ax.bar(labels, counts)
    ax.set_ylim(0, max(counts) * 1.15)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    height + 1,
                    str(count),
                    ha='center',
                    va='bottom',
                    fontsize=8)
    ax.set_title("E-Value")
    ax.set_xlabel("E-Value range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()
    
    
    
def draw_alength_histogram(chart_widget, lengths, bin_count=5):
    if len(lengths) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Alignment Lengths Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_len = min(lengths)
    max_len = max(lengths)
    bin_size = max(1, round((max_len - min_len) / bin_count))

    bins = list(range(min_len, max_len + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)
    
    

    for length in lengths:
        for i in range(len(bins) - 1):
            if bins[i] <= length < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if length == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bar_container = ax.bar(labels, counts, picker=True)
    ax.set_ylim(0, max(counts) * 1.15)
    for rect, count in zip(bar_container, counts):
        height = rect.get_height()
        if height > 0:
            ax.text(rect.get_x() + rect.get_width() / 2, height + 1,
                    str(count), ha='center', va='bottom', fontsize=8)
    ax.set_title("Alignment Length")
    ax.set_xlabel("Length range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()

    parent = chart_widget.parentWidget()
    if parent:
        parent.length_bins = bins
        parent.length_bar_patches = list(bar_container.patches)
        
        
        
def draw_identities_histogram(chart_widget, headers):
    import re

    scores = []
    for header in headers:
        match = re.search(r"Identities: ([\d\.]+)", header)
        if match:
            try:
                scores.append(float(match.group(1)))
            except ValueError:
                pass

    if not scores or len(scores) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Identity Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_score = min(scores)
    max_score = max(scores)
    bin_count = 5
    bin_size = max(1, round((max_score - min_score) / bin_count))
    bins = list(range(int(min_score), int(max_score) + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for score in scores:
        for i in range(len(bins) - 1):
            if bins[i] <= score < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if score == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bars = ax.bar(labels, counts)
    ax.set_ylim(0, max(counts) * 1.15)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    height + 1,
                    str(count),
                    ha='center',
                    va='bottom',
                    fontsize=8)
    
                    
    ax.set_title("Identity")
    ax.set_xlabel("Ident. range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()
    
    

def draw_positives_histogram(chart_widget, headers):
    import re

    scores = []
    for header in headers:
        match = re.search(r"Positives: ([\d\.]+)", header)
        if match:
            try:
                scores.append(float(match.group(1)))
            except ValueError:
                pass

    if not scores or len(scores) < 2:
        chart_widget.figure.clear()
        ax = chart_widget.figure.add_subplot(111)
        ax.set_title("No Similarity Info", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_score = min(scores)
    max_score = max(scores)
    bin_count = 5
    bin_size = max(1, round((max_score - min_score) / bin_count))
    bins = list(range(int(min_score), int(max_score) + bin_size, bin_size))
    counts = [0] * (len(bins) - 1)

    for score in scores:
        for i in range(len(bins) - 1):
            if bins[i] <= score < bins[i + 1]:
                counts[i] += 1
                break
        else:
            if score == bins[-1]:
                counts[-1] += 1

    chart_widget.figure.clear()
    ax = chart_widget.figure.add_subplot(111)
    labels = [f"{round(bins[i])}-{round(bins[i+1]-1)}" for i in range(len(bins) - 1)]
    bars = ax.bar(labels, counts)
    ax.set_ylim(0, max(counts) * 1.15)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    height + 1,
                    str(count),
                    ha='center',
                    va='bottom',
                    fontsize=8)
    
                    
    ax.set_title("Similarity")
    ax.set_xlabel("Sim. range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    fig = chart_widget.figure
    fig.tight_layout()
    chart_widget.canvas.draw()
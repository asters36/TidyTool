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
    ax.set_title("Lengths Histogram")
    ax.set_xlabel("Length range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
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
        ax.set_title("NO DATA")
        ax.set_xticks([])  # brak osi X
        ax.set_yticks([])  # brak osi Y
        ax.axis('off')
        chart_widget.canvas.draw()
        return

    min_score = min(scores)
    max_score = max(scores)
    bin_count = 4
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
    ax.bar(labels, counts)
    ax.set_title("Bit Score Histogram")
    ax.set_xlabel("Score range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
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
        ax.set_title("NO DATA")
        ax.set_xticks([])  # brak osi X
        ax.set_yticks([])  # brak osi Y
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
    ax.bar(labels, counts)
    ax.set_title("E-Values Histogram")
    ax.set_xlabel("E-Value range")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=0, labelsize=8)
    chart_widget.canvas.draw()
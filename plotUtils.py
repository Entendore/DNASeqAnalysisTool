import matplotlib.pyplot as plt
import seaborn as sns
from analysis import calculate_base_frequencies
import plotly.graph_objects as go

def plot_base_distribution(freqs, output_file=None, svg=False):
    if svg and output_file:
        sns.barplot(x=list(freqs.keys()), y=list(freqs.values()))
        plt.title("Base Frequency Distribution")
        plt.savefig(output_file)
        plt.close()
    else:
        fig = go.Figure([go.Bar(x=list(freqs.keys()), y=list(freqs.values()))])
        fig.update_layout(title="Base Frequency Distribution", xaxis_title="Base", yaxis_title="Frequency")
        if output_file:
            fig.write_html(output_file.replace('.svg', '.html').replace('.png', '.html'))
        else:
            fig.show()

def plot_gc_profile(seq, output_file=None, window_size=100, svg=False):
    gc_values = []
    for i in range(0, len(seq) - window_size + 1, window_size):
        window = seq[i:i + window_size]
        _, gc = calculate_base_frequencies(window)
        gc_values.append(gc)
    if svg and output_file:
        plt.plot(gc_values)
        plt.title("GC Content Profile")
        plt.xlabel("Window")
        plt.ylabel("GC Content")
        plt.savefig(output_file)
        plt.close()
    else:
        fig = go.Figure(go.Scatter(y=gc_values))
        fig.update_layout(title="GC Content Profile", xaxis_title="Window", yaxis_title="GC Content")
        if output_file:
            fig.write_html(output_file.replace('.svg', '.html').replace('.png', '.html'))
        else:
            fig.show()

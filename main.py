import os
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from analysis import calculate_base_frequencies, find_orfs
from utils import load_sequences
from plotUtils import plot_base_distribution, plot_gc_profile

os.makedirs("outputs", exist_ok=True)
os.makedirs("plots", exist_ok=True)

# Change here to load your input file path and format
input_file = "data/example.fasta"
input_format = "fasta"

seq_records = load_sequences(input_file, file_format=input_format)

all_orfs = []
all_freqs = {}

for record in seq_records:
    freqs, gc = calculate_base_frequencies(record.seq)
    all_freqs[record.id] = (freqs, gc)
    orfs = find_orfs(record, min_length=100, table=1)
    all_orfs.extend(orfs)

    plot_base_distribution(freqs, f"plots/base_distribution_{record.id}.svg", svg=True)
    plot_gc_profile(record.seq, f"plots/gc_profile_{record.id}.svg", svg=True)

df = pd.DataFrame(all_orfs, columns=["SequenceID", "Strand", "Start", "End", "Length"])
df.to_csv("outputs/orf_report.csv", index=False)

# Export ORFs as FASTA
output_records = []
for sid, strand, start, end, _ in all_orfs:
    rec = next(r for r in seq_records if r.id == sid)
    if strand == 1:
        orf_seq = rec.seq[start:end]
    else:
        orf_seq = rec.seq[start:end].reverse_complement()
    output_records.append(SeqRecord(orf_seq, id=f"{sid}_orf_{start}_{end}", description=""))
SeqIO.write(output_records, "outputs/orfs.fasta", "fasta")

print(f"Processed {len(seq_records)} sequences")
print(f"Total ORFs found: {len(all_orfs)}")
print("Reports and plots saved in outputs/ and plots/")

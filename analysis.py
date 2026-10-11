from collections import Counter
from Bio.Seq import Seq

def calculate_base_frequencies(seq):
    seq = str(seq).upper()
    counts = Counter(seq)
    total = sum(counts[base] for base in "ATGC")
    freqs = {base: counts.get(base, 0) / total for base in "ATGC"}
    gc_content = (counts.get("G", 0) + counts.get("C", 0)) / total
    return freqs, gc_content

def find_orfs(seq_record, min_length=100, table=1):
    """
    Find ORFs in a sequence record using given genetic code table.
    Returns list of (seq_id, strand, start, end, length).
    """
    orfs = []
    seq = seq_record.seq
    for strand, nucleotide_seq in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = nucleotide_seq[frame:].translate(table=table, to_stop=False)
            aa_seq = str(trans)
            start = 0
            while start < len(aa_seq):
                start_codon = aa_seq.find("M", start)
                if start_codon == -1:
                    break
                stop_codon = aa_seq.find("*", start_codon)
                if stop_codon == -1:
                    break
                length = stop_codon - start_codon
                if length >= min_length:
                    orf_start = frame + start_codon * 3
                    orf_end = frame + stop_codon * 3 + 3
                    if strand == -1:
                        orf_start, orf_end = len(seq) - orf_end, len(seq) - orf_start
                    orfs.append((seq_record.id, strand, orf_start, orf_end, length))
                start = stop_codon + 1
    return orfs

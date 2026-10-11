from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def load_fasta(filepath):
    if isinstance(filepath, str):
        with open(filepath, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
    else:
        # Assume file-like object (e.g., Streamlit upload)
        records = list(SeqIO.parse(filepath, "fasta"))
    return records

def save_fasta(records, filepath):
    with open(filepath, "w") as f:
        SeqIO.write(records, f, "fasta")

def load_sequences(file, file_format="fasta"):
    """
    Load sequences from a file or file-like object, supporting multiple formats.
    """
    file_format = file_format.lower()
    if file_format in ["fasta", "genbank"]:
        return list(SeqIO.parse(file, file_format))
    elif file_format == "raw_text":
        text = file.read().decode("utf-8").strip().replace("\n", "").upper()
        return [SeqRecord(Seq(text), id="raw_seq")]
    elif file_format == "csv":
        # file is file-like object, read with pandas
        df = pd.read_csv(file)
        # Expect user to handle selecting column in app.py
        return df
    else:
        raise ValueError(f"Unsupported file format: {file_format}")

def save_sequences(sequences, file, file_format="fasta"):
    with open(file, "w") as f:
        SeqIO.write(sequences, f, file_format)

def load_raw_text(file):
    text = file.read().decode("utf-8").strip().replace("\n", "").upper()
    return [SeqRecord(Seq(text), id="raw_seq")]

def save_raw_text(text, file):
    with open(file, "w") as f:
        f.write(text)


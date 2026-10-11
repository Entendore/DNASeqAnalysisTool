import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez, SeqIO
from Bio.Seq import translate
from analysis import calculate_base_frequencies, find_orfs
from utils import load_sequences
import os
import plotly.express as px
import io
import base64

st.title("DNA Sequence Analysis Tool")

os.makedirs("outputs", exist_ok=True)

def get_download_link(df, filename):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download {filename}</a>'
    return href

def get_fasta_download_link(records, filename):
    fasta_io = io.StringIO()
    SeqIO.write(records, fasta_io, "fasta")
    fasta_str = fasta_io.getvalue()
    b64 = base64.b64encode(fasta_str.encode()).decode()
    href = f'<a href="data:file/fasta;base64,{b64}" download="{filename}">Download {filename}</a>'
    return href

format_choice = st.selectbox("Select file format", options=["fasta", "genbank", "raw_text", "csv", "ncbi_accession"])

sequences = []

try:
    if format_choice == "ncbi_accession":
        Entrez.email = st.text_input("Enter your email (required for NCBI fetch)")
        accession = st.text_input("Enter NCBI accession number")
        if st.button("Fetch from NCBI"):
            if Entrez.email and accession:
                with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
                    sequences = list(SeqIO.parse(handle, "fasta"))
                st.success(f"Fetched {len(sequences)} sequence(s)")
            else:
                st.error("Please enter a valid email and accession number.")

    elif format_choice == "raw_text":
        raw_seq = st.text_area("Paste raw nucleotide sequence here")
        if raw_seq:
            sequences = [SeqRecord(Seq(raw_seq.strip().replace("\n", "").upper()), id="raw_seq")]

    elif format_choice == "csv":
        csv_file = st.file_uploader("Upload CSV file with sequences", type=["csv"])
        if csv_file:
            df = pd.read_csv(csv_file)
            seq_col = st.selectbox("Select column with sequences", df.columns)
            if seq_col:
                sequences = [SeqRecord(Seq(seq.upper()), id=str(i)) for i, seq in enumerate(df[seq_col])]

    else:
        uploaded_file = st.file_uploader(f"Upload {format_choice} file", type=[format_choice, format_choice + "k"])
        if uploaded_file:
            sequences = load_sequences(uploaded_file, file_format=format_choice)

except Exception as e:
    st.error(f"Error loading sequences: {e}")

window_size = st.slider("GC Content Window Size", min_value=10, max_value=500, value=100, step=10)
min_orf_length = st.number_input("Minimum ORF Length (AA)", min_value=10, max_value=1000, value=100, step=10)
translation_table = st.selectbox("Translation Table (NCBI)", options=[1, 2, 3, 4, 5, 11, 12], index=0,
                                 help="Select NCBI translation table number")

if sequences:
    st.write(f"Processing {len(sequences)} sequence(s)...")

    all_orfs = []
    orf_details = []

    for record in sequences:
        st.subheader(f"Sequence ID: {record.id}")

        try:
            freqs, gc = calculate_base_frequencies(record.seq)
            st.write(f"Base Frequencies: {freqs}")
            st.write(f"GC Content: {gc:.2%}")

            orfs = find_orfs(record, min_length=min_orf_length, table=translation_table)
            all_orfs.extend(orfs)

            # Interactive base freq plot
            fig1 = px.bar(x=list(freqs.keys()), y=list(freqs.values()),
                          labels={"x": "Base", "y": "Frequency"},
                          title=f"Base Frequency Distribution: {record.id}")
            st.plotly_chart(fig1)

            # GC content profile plot
            gc_values = []
            for i in range(0, len(record.seq) - window_size + 1, window_size):
                window = record.seq[i:i + window_size]
                _, gc_w = calculate_base_frequencies(window)
                gc_values.append(gc_w)
            fig2 = px.line(y=gc_values, labels={"y": "GC Content", "x": "Window"},
                           title=f"GC Content Profile: {record.id}")
            st.plotly_chart(fig2)

            # Build detailed ORF info with AA sequences
            for sid, strand, start, end, length in orfs:
                if strand == 1:
                    orf_seq = record.seq[start:end]
                else:
                    orf_seq = record.seq[start:end].reverse_complement()
                aa_seq = orf_seq.translate(table=translation_table)
                orf_details.append({
                    "SequenceID": sid,
                    "Strand": strand,
                    "Start": start,
                    "End": end,
                    "Length": length,
                    "AA_Sequence": str(aa_seq)
                })

        except Exception as e:
            st.error(f"Error processing sequence {record.id}: {e}")

    st.write(f"Total ORFs found: {len(all_orfs)}")

    if orf_details:
        df = pd.DataFrame(orf_details)
        st.dataframe(df)

        # Download links
        st.markdown(get_download_link(df, "orf_report.csv"), unsafe_allow_html=True)

        output_records = []
        for item in orf_details:
            seq_id = item["SequenceID"]
            strand = item["Strand"]
            start = item["Start"]
            end = item["End"]
            rec = next(r for r in sequences if r.id == seq_id)
            if strand == 1:
                orf_seq = rec.seq[start:end]
            else:
                orf_seq = rec.seq[start:end].reverse_complement()
            output_records.append(SeqRecord(orf_seq, id=f"{seq_id}_orf_{start}_{end}", description=""))

        st.markdown(get_fasta_download_link(output_records, "orfs.fasta"), unsafe_allow_html=True)

else:
    st.info("Upload or fetch sequence data to get started.")

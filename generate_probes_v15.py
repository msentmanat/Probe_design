import sys
import pandas as pd
import re
import subprocess
import tempfile
import os
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
from snapgene_reader import snapgene_file_to_dict

def has_homopolymer(primer):
    return bool(re.search(r"(A{4,}|T{4,}|C{4,}|G{4,})", primer.upper()))

def ends_in_gc(primer):
    return primer.upper()[-1] in "GC"

def is_valid_primer(primer):
    primer = primer.upper()
    return not has_homopolymer(primer) and ends_in_gc(primer)

def get_feature_for_position(feature_list, pos):
    for feat in feature_list:
        f_start, f_end = feat['range']
        if f_start <= pos <= f_end:
            return feat['name']
    return None

def is_within_feature_range(feature_list, start_pos, end_pos):
    for feat in feature_list:
        f_start, f_end = feat['range']
        if f_start <= start_pos and end_pos <= f_end:
            return feat['name']
    return None

def is_position_in_feature(feature_list, pos):
    for feat in feature_list:
        f_start, f_end = feat['range']
        if f_start <= pos <= f_end:
            return feat['name']
    return None

def percent_masked(seq):
    return (sum(1 for c in seq if c.islower()) / len(seq)) * 100
    
def mask_ns(sequence):
    """Mask ambiguous bases ('N') in a sequence by converting them to lowercase."""
    return ''.join([base.lower() if base.upper() == 'N' else base.upper() for base in sequence])


def run_repeatmasker(sequence, tempdir, species):
    fasta_path = os.path.join(tempdir, "input.fa")
    with open(fasta_path, "w") as f:
        f.write(">sequence\n")
        f.write(sequence + "\n")

    try:
        subprocess.run(["RepeatMasker", "-nolow", "-no_is", "-norna", "-xsmall", "-species", species, fasta_path],
               check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
               
    except subprocess.CalledProcessError as e:
        print("❌ RepeatMasker failed:", e.stderr.decode())
        sys.exit(1)

    masked_path = os.path.join(tempdir, "input.fa.masked")
    if not os.path.exists(masked_path):
        print("❌ Masked output not found.")
        sys.exit(1)

    masked_record = next(SeqIO.parse(masked_path, "fasta"))
    return str(masked_record.seq)
    
def run_repeatmasker_with_engine(sequence, tempdir, species, engine):
    fasta_path = os.path.join(tempdir, f"input_{engine}.fa")
    with open(fasta_path, "w") as f:
        f.write(">sequence\n")
        f.write(sequence + "\n")

    try:
        subprocess.run(["RepeatMasker", "-engine", engine, "-nolow", "-no_is", "-norna",
                        "-xsmall", "-species", species, fasta_path],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"❌ RepeatMasker ({engine}) failed:", e.stderr.decode())
        sys.exit(1)

    masked_path = f"{fasta_path}.masked"
    if not os.path.exists(masked_path):
        print(f"❌ Masked output ({engine}) not found.")
        sys.exit(1)

    masked_record = next(SeqIO.parse(masked_path, "fasta"))
    return str(masked_record.seq)

def merge_maskings(original, masked1, masked2):
    merged = []
    for o, m1, m2 in zip(original, masked1, masked2):
        if m1.islower() or m2.islower():
            merged.append(o.lower())
        else:
            merged.append(o.upper())
    return ''.join(merged)




def find_primers(sequence, masked_sequence, features, project_name,
                 min_amp=126, max_amp=326, ideal_amp=226,
                 min_primer_len=20, max_primer_len=24, min_spacing=650,
                 min_gc=30, max_gc=80, gc_tolerance=5,
                 max_masked_percent=100):

    primer_sets = []
    seq_len = len(sequence)
    index = 1
    start = 0
    forward_tag = "ggcatttcagtcaggtgcccaatgtacc"
    reverse_tag = "gttccgggtaggcagttcgctccaagct"

    while start + min_amp <= seq_len:
        forward_primer = None
        f_gc = None
        f_len = None

        for fwd_len in range(min_primer_len, max_primer_len + 1):
            if start + fwd_len > seq_len:
                break

            candidate_fwd = sequence[start:start + fwd_len].upper()
            candidate_masked_fwd = masked_sequence[start:start + fwd_len]
            if percent_masked(candidate_masked_fwd) > max_masked_percent:
                continue

            candidate_f_gc = gc_fraction(candidate_fwd) * 100

            if (min_gc <= candidate_f_gc <= max_gc) and is_valid_primer(candidate_fwd):
                forward_primer = candidate_fwd
                f_gc = candidate_f_gc
                f_len = fwd_len
                break

        if forward_primer is None:
            start += 1
            continue

        search_positions = []
        offset = 0
        while True:
            down = ideal_amp - offset
            up = ideal_amp + offset
            if min_amp <= down <= max_amp:
                search_positions.append(down)
            if min_amp <= up <= max_amp and up != down:
                search_positions.append(up)
            if down <= min_amp and up >= max_amp:
                break
            offset += 1

        reverse_primer = None
        best_r_gc = None
        best_end = None
        best_r_len = None

        for amp_size in search_positions:
            search_start = start + amp_size
            if search_start >= seq_len:
                continue

            for r_len in range(min_primer_len, max_primer_len + 1):
                end = search_start + r_len
                if end > seq_len:
                    continue

                rev_start = end - r_len
                candidate_rev = Seq(sequence[rev_start:end]).reverse_complement()
                candidate_rev_str = str(candidate_rev).upper()

                masked_candidate_rev = masked_sequence[rev_start:end]
                if percent_masked(masked_candidate_rev) > max_masked_percent:
                    continue

                candidate_r_gc = gc_fraction(candidate_rev_str) * 100

                if (min_gc <= candidate_r_gc <= max_gc and
                    abs(f_gc - candidate_r_gc) <= gc_tolerance and
                    is_valid_primer(candidate_rev_str)):

                    reverse_primer = candidate_rev_str
                    best_r_gc = candidate_r_gc
                    best_end = end
                    best_r_len = r_len
                    break
            if reverse_primer:
                break

        if reverse_primer is None:
            start += 1
            continue

        amplicon = sequence[start:best_end].upper()
        masked_amplicon = masked_sequence[start:best_end]
        if percent_masked(masked_amplicon) > max_masked_percent:
            start += 1
            continue

        amp_gc = gc_fraction(amplicon) * 100

        fwd_start = start + 1
        fwd_end = fwd_start + f_len - 1
        rev_end = best_end
        rev_start = rev_end - best_r_len + 1

        fwd_feature = is_within_feature_range(features, fwd_start, fwd_end)
        rev_feature = is_position_in_feature(features, rev_end)

        if fwd_feature == rev_feature and fwd_feature is not None:
            feature_name = fwd_feature
        else:
            feature_name = ""
        
        suffix = ""
        #suffix = f".{feature_name}" if feature_name else ""

        primer_sets.append({
            "Primer Name": f"{project_name}{suffix}.F{index}",
            "Forward Primer": forward_primer,
            "Reverse Primer": "",
            "Primer GC%": f_gc,
            "Primer Length": f_len,
            "Amplicon Size": len(amplicon),
            "Amplicon GC%": amp_gc,
            "Primer to Order": forward_tag + forward_primer,
            "Flag": forward_tag
        })

        primer_sets.append({
            "Primer Name": f"{project_name}{suffix}.R{index}",
            "Forward Primer": "",
            "Reverse Primer": reverse_primer,
            "Primer GC%": best_r_gc,
            "Primer Length": best_r_len,
            "Amplicon Size": "",
            "Amplicon GC%": "",
            "Primer to Order": reverse_tag + reverse_primer,
            "Flag": reverse_tag
        })

        index += 1
        start = best_end + min_spacing

    return primer_sets

if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser(description="Design probes from a SnapGene file with optional RepeatMasker masking.")
    parser.add_argument("project_name", help="Project name")
    parser.add_argument("input_file", help="Input .dna SnapGene file")
    parser.add_argument("--use_repeatmasker", action="store_true", help="Use RepeatMasker to mask repeats")
    parser.add_argument("--species", default="mouse", help="Species name for RepeatMasker (e.g., mouse, human, etc.)")
    parser.add_argument("--mask_ns", action="store_true", help="Mask ambiguous bases (N) as repeats")

    
    args = parser.parse_args()
    
    project_name = args.project_name
    snapgene_file = args.input_file
    use_repeatmasker = args.use_repeatmasker
    species = args.species

    snap = snapgene_file_to_dict(snapgene_file)
    sequence = snap['seq'].upper()


    if use_repeatmasker:
        with tempfile.TemporaryDirectory() as tempdir:
            masked_rmblast = run_repeatmasker_with_engine(sequence, tempdir, species, engine="rmblast")
            #masked_hmmer = run_repeatmasker_with_engine(sequence, tempdir, species, engine="hmmer")
            #masked_sequence = merge_maskings(sequence, masked_rmblast, masked_hmmer)
            masked_sequence = masked_rmblast
            masked_fasta_path = f"{project_name}_masked_sequence.fa"
            with open(masked_fasta_path, "w") as masked_out:
                masked_out.write(f">{project_name}_masked\n{masked_sequence}\n")
            print(f"✅ Combined masked sequence saved to {masked_fasta_path}")


    else:
        print("⚠️ Skipping RepeatMasker. Using unmasked sequence.")
        masked_sequence = sequence
    
    if args.mask_ns:
        masked_sequence = mask_ns(masked_sequence)
        print("✅ Applied N-base masking to sequence.")



    # Extract features, skip "source"
    raw_features = snap.get('features', [])
    features = []
    sequence_length = len(sequence)

    for feat in raw_features:
        feat_type = feat.get("type", "").lower()
        name = feat.get("name", "Unnamed")

        # Skip full-sequence "source" features (likely metadata)
        skip = False
        for seg in feat.get("segments", []):
            range_str = seg.get("@range")
            if range_str:
                match = re.match(r"(\d+)-(\d+)", range_str)
                if match:
                    start, end = int(match.group(1)), int(match.group(2))
                    if start == 1 and end == sequence_length:
                        skip = True
        if skip:
            continue
    
        # Otherwise, extract segments
        for seg in feat.get("segments", []):
            range_str = seg.get("@range")
            if range_str:
                match = re.match(r"(\d+)-(\d+)", range_str)
                if match:
                    start, end = int(match.group(1)), int(match.group(2))
                    features.append({
                        "name": name,
                        "range": (start, end)
                    })



    print("Parsed features:")
    for f in features:
        print(f" - {f['name']} {f['range']}")




    print(f"✅ Parsed {len(features)} annotated features.")
    primer_data = find_primers(sequence, masked_sequence, features, project_name)

    file_name = f"{project_name}_probes.csv"
    file_name2 = f"{project_name}_primers.csv"
    df = pd.DataFrame(primer_data)
    df["Placeholder"] = ""
    column_order = ["Primer Name", "Primer to Order", "Placeholder", "Amplicon Size", "Amplicon GC%", "Flag",
                    "Forward Primer", "Reverse Primer", "Primer GC%", "Primer Length"]
    df = df[column_order]
    df.to_csv(file_name, index=False)

    print(f"\n✅ Probes saved to {file_name}")


#    # -------------------------------------------------------
#    # Append primer-pair table to the same probes CSV
#    # -------------------------------------------------------
#
#    rows = list(df.iterrows())

#    with open(file_name, "a") as f_out:
#        f_out.write("\n")
#        f_out.write("Primer pair,Forward primer,Reverse primer\n")
#
#        for (_, r1), (_, r2) in zip(rows[::2], rows[1::2]):
#
#            name_f = r1["Primer Name"]
#            name_r = r2["Primer Name"]
#
#            # Build pair name like Project.F1+R1
#            base = name_f.rsplit(".F", 1)[0]
#            f_index = name_f.rsplit(".F", 1)[1]
#            r_index = name_r.rsplit(".R", 1)[1]
#
#            pair_name = f"{base}.F{f_index}+R{r_index}"
#
#            f_out.write(
#                f"{pair_name},{r1['Primer to Order']},{r2['Primer to Order']}\n"
#            )


# Append a blank line and simplified Primer Name + Sequence table
    with open(file_name, "a") as f_out:
        f_out.write("\n")
        f_out.write("Primer Name,Primer Sequence\n")
    
        for _, row in df.iterrows():
            name = row["Primer Name"]
            seq = row["Forward Primer"] or row["Reverse Primer"]
            if seq:
                f_out.write(f"{name},{seq}\n")




#    with open(file_name2, "a") as f_out2:
#        for _, row in df.iterrows():
#            name = row["Primer Name"]
#            seq = row["Forward Primer"] or row["Reverse Primer"]
#            if seq:
#                f_out2.write(f"{name},{seq}\n")


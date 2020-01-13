from Bio import SeqIO
from numpy import diff, where, split
import argparse
import pandas as pd
from wiggle_parser import WiggleParser as wp
import glob

def main():
    # Params
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    parser.add_argument("--pre_signal_offset", required=True, help="", type=int)
    parser.add_argument("--post_signal_offset", required=True, help="", type=int)
    parser.add_argument("--min_coverage", required=True, help="", type=float)
    parser.add_argument("--max_interruption", required=True, help="", type=int)
    parser.add_argument("--window_size", required=True, help="", type=int)
    parser.add_argument("--tolerance", required=True, help="", type=int)
    parser.add_argument("--merge_range", default=0, required=False, help="", type=int)
    parser.add_argument("--min_len", default=None, required=False, help="", type=int)
    parser.add_argument("--base", default="T", required=False, help="", type=str)
    args = parser.parse_args()
    # ---------------------------
    if args.min_len <= args.window_size or args.min_len is None:
        print("WARNING: Minimum length will be overridden by the window size.")
        args.min_len = args.window_size
    print("Loading sequence file...")
    fasta_parsed = SeqIO.parse(glob.glob(args.fasta_in)[0], "fasta")
    wig_files = glob.glob(args.wigs_in)
    f_wigs_parsed, r_wigs_parsed = merge_wigs_by_max(wig_files, fasta_parsed)
    accession = ""
    ret_list = []
    counters = {}
    # The following line is repeated due to list exhaustion
    fasta_parsed = SeqIO.parse(glob.glob(args.fasta_in)[0], "fasta")
    for seq_record in fasta_parsed:
        f_seq_str = str(seq_record.seq)
        print(f_seq_str)
        accession = seq_record.id
        f_positions, r_positions = group_positions(f_seq_str, args.base, args.max_interruption, args.window_size,
                                                   args.tolerance, args.min_len)
        # f_positions, r_positions = seek_window(f_seq_str, args.window_size, args.tolerance)
        counters[f'wig_pos_count_{accession}'] = 0
        f_wig_df_sliced = f_wigs_parsed[accession][f_wigs_parsed[accession][1] >= args.min_coverage]
        r_wig_df_sliced = r_wigs_parsed[accession][r_wigs_parsed[accession][1] <= args.min_coverage*-1]
        counters[f'wig_pos_count_{accession}'] += f_wig_df_sliced.shape[0] + r_wig_df_sliced.shape[0]

        for i in f_positions:
            i[0] -= args.pre_signal_offset
            i[1] += args.post_signal_offset
        for i in r_positions:
            i[0] -= args.post_signal_offset
            i[1] += args.pre_signal_offset
        # Merge overlapping position after adding the offsets within a range
        f_positions = merge_interval_lists(f_positions, args.merge_range)
        r_positions = merge_interval_lists(r_positions, args.merge_range)

        counters[f'poly_signals_count_{accession}'] = len(f_positions) + len(r_positions)
        counters[f'a_cov_drop_matches_{accession}'] = 0
        counters[f'f_cov_drop_matches_{accession}'] = 0
        counters[f'r_cov_drop_matches_{accession}'] = 0
        counters[f'pos_not_in_cov_{accession}'] = 0
        counters[f'pos_no_match_{accession}'] = 0
        for position in f_positions:
            if get_score_of_wig_loc(f_wig_df_sliced, list(range(position[0], position[1]))):
                counters[f'f_cov_drop_matches_{accession}'] += 1
                ret_list.append([accession, position[0], position[1], "+"])
            else:
                if get_score_of_wig_loc(f_wigs_parsed[accession], list(range(position[0], position[1]))):
                    counters[f'pos_no_match_{accession}'] += 1
                else:
                    counters[f'pos_not_in_cov_{accession}'] += 1
        for position in r_positions:
            if get_score_of_wig_loc(r_wig_df_sliced, list(range(position[0], position[1]))):
                counters[f'r_cov_drop_matches_{accession}'] += 1
                ret_list.append([accession, position[0], position[1], "-"])
            else:
                if get_score_of_wig_loc(r_wigs_parsed[accession], list(range(position[0], position[1]))):
                    counters[f'pos_no_match_{accession}'] += 1
                else:
                    counters[f'pos_not_in_cov_{accession}'] += 1
        counters[f'a_cov_drop_matches_{accession}'] = counters[f'f_cov_drop_matches_{accession}'] + \
                                                      counters[f'r_cov_drop_matches_{accession}']
    # OUTPUT
    print(f"Given parameters:\n"
          f"\t- Pre-signal offset: {args.pre_signal_offset}\n"
          f"\t- Post-signal offset: {args.post_signal_offset}\n"
          f"\t- Minimum coverage: {args.min_coverage}\n"
          f"\t- Maximum interruption: {args.max_interruption}\n"
          f"\t- Window size: {args.window_size}\n"
          f"\t- Tolerance: {args.tolerance}\n"
          f"\t- Merge range: {args.merge_range}\n"
          f"\t- Minimum poly-{args.base} length: {args.min_len}\n"
          f"\t- Base: {args.base}")
    print(f"Output:\n"
          f"\t- Total coverage drop matches (sum): "
          f"{sum(v for k, v in counters.items() if 'a_cov_drop_matches_' in k):,}\n"
          f"\t- Total coverage drop matches in forward: "
          f"{sum(v for k, v in counters.items() if 'f_cov_drop_matches_' in k):,}\n"
          f"\t- Total coverage drop matches in reverse: "
          f"{sum(v for k, v in counters.items() if 'r_cov_drop_matches_' in k):,}\n"
          f"\t- Total Poly-{args.base} signals: {sum(v for k, v in counters.items() if 'poly_signals_count_' in k):,}")
    print(f"\t- Poly-{args.base} signals that has no coverage: "
          f"{sum(v for k, v in counters.items() if 'pos_no_match_' in k):,}\n"
          f"\t- Poly-{args.base} signals that found in the coverage but does not meet the matching parameters: "
          f"{sum(v for k, v in counters.items() if 'pos_not_in_cov_' in k):,}\n"
          f"\t- Total number of positions in coverage: "
          f"{sum(v for k, v in counters.items() if 'wig_pos_count_' in k):,}")

    print("Writing GFF file...")
    term_gff_str = ""
    count = 0
    current_accession = ""
    out_df = pd.DataFrame.from_records(ret_list)
    out_df = out_df.sort_values([0, 1])
    for index, row in out_df.iterrows():
        if current_accession != row[0] or current_accession == "":
            # term_gff_str += "###\n"
            current_accession = row[0]
            count = 0
        count += 1
        term_gff_str += \
            f"{row[0]}\t" + \
            f"term_last_base_finder\t" + \
            f"possible_term_end\t" + \
            f"{row[1]}\t" + \
            f"{row[2]}\t" + \
            f".\t" + \
            f"{row[3]}\t" + \
            f".\t" + \
            f"id={current_accession}_{row[3]}_term_end{count};" + \
            f"name={current_accession}_{row[3]}_term_end{count}\n"
    outfile = open(args.gff_out, "w")
    outfile.write(f"###gff-version 3\n{term_gff_str}###")
    outfile.close()


def merge_wigs_by_max(wig_files, fasta_parsed):
    # Handling multiple wig files
    f_wigs_parsed = {}
    r_wigs_parsed = {}
    for seq_record in fasta_parsed:
        f_wigs_parsed[seq_record.id] = pd.DataFrame(data=range(1, len(seq_record.seq), 1))
        r_wigs_parsed[seq_record.id] = pd.DataFrame(data=range(1, len(seq_record.seq), 1))

    for wig in wig_files:
        x = wp(wig).parse()
        for key, value in x.items():
            if key in r_wigs_parsed.keys():
                if value[value[1] < 0].empty:
                    f_wigs_parsed[key] = pd.merge(how='outer', left=f_wigs_parsed[key], right=value, left_on=0,
                                                  right_on=0).fillna(0.0)
                if value[value[1] > 0].empty:
                    r_wigs_parsed[key] = pd.merge(how='outer', left=r_wigs_parsed[key], right=value, left_on=0,
                                                  right_on=0).fillna(0.0)
    for accession in f_wigs_parsed.keys():
        f_wigs_parsed[accession][1] = f_wigs_parsed[accession].iloc[:, 1:-1].max(axis=1)
        f_wigs_parsed[accession] = f_wigs_parsed[accession].iloc[:, [0, -1]]
    for accession in r_wigs_parsed.keys():
        r_wigs_parsed[accession][1] = r_wigs_parsed[accession].iloc[:, 1:-1].min(axis=1)
        r_wigs_parsed[accession] = r_wigs_parsed[accession].iloc[:, [0, -1]]
    return f_wigs_parsed, r_wigs_parsed


def drop_invalid_signals(all_signals, window_size, tolerance):
    clean_signals = []
    for signal in all_signals:
        # Drop any signal shorter than the window
        if len(signal) < window_size - 1 - tolerance:
            continue
        # Seek the window across the signal
        for index, pos in enumerate(signal):
            # first check if the end of window does not exceed the list size
            if len(signal) < index + window_size:
                break
            # Check if the sliding window contains the required positions
            if 0 <= (signal[index + (window_size - 1)] - pos) - window_size <= tolerance:
                clean_signals.append([signal[0], signal[-1]])
                break
    return clean_signals


def group_positions(seq_str, base, max_interruption, window_size, tolerance, min_len=None):
    complement_base = lambda x: "T" if base == "A" else "A"
    f_indices = [i for i, a in enumerate(seq_str, 1) if a == base]
    r_indices = [i for i, a in enumerate(seq_str, 1) if a == complement_base(base)]
    f_indices.sort()
    r_indices.sort()
    # Get all signals to any length with max interruption
    f_signals = list(map(list, split(f_indices, where(diff(f_indices) > max_interruption)[0] + 1)))
    r_signals = list(map(list, split(r_indices, where(diff(r_indices) > max_interruption)[0] + 1)))
    f_signals = drop_invalid_signals(f_signals, window_size, tolerance)
    r_signals = drop_invalid_signals(r_signals, window_size, tolerance)
    f_poly_base_signal_locations = [[i[0], i[-1]] for i in f_signals if i[-1] - i[0] >= min_len]
    r_poly_base_signal_locations = [[i[0], i[-1]] for i in r_signals if i[-1] - i[0] >= min_len]
    return f_poly_base_signal_locations, r_poly_base_signal_locations


def merge_interval_lists(list_in, merge_range):
    list_out = []
    for loc in list_in:
        if len(list_out) == 0:
            list_out.append(loc)
        else:
            if loc[0] in range(list_out[-1][0], list_out[-1][-1] + merge_range):
                list_out[-1][-1] = loc[-1]
            else:
                list_out.append(loc)
    return list_out


def seek_window(seq_str, window_size, tolerance):
    f_locations = []
    r_locations = []
    for index, nt in enumerate(seq_str):
        if len(seq_str) < index + window_size:
            break
        if 0 <= window_size - seq_str[index:index + window_size].count("T") <= tolerance:
            f_locations.append([index + 1, index + 1 + window_size])
        if 0 <= window_size - seq_str[index:index + window_size].count("A") <= tolerance:
            r_locations.append([index + 1, index + 1 + window_size])
    f_locations = merge_interval_lists(f_locations, 0)
    r_locations = merge_interval_lists(r_locations, 0)
    return f_locations, r_locations


def get_score_of_wig_loc(wig_df, pos):
    x = wig_df[wig_df[0].isin(range(pos[0], pos[-1]))]
    if x.empty:
        return False
    else:
        return True


main()

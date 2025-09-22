import pandas as pd
from collections import defaultdict
from math import isclose
import numpy as np

def mutation_info(seq1, seq2, hotspot_sites):
    total_mut = 0
    hotspot_mut = 0
    adjusted_hotspots = set([pos - 1 for pos in hotspot_sites])  # convert 1-based to 0-based
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            total_mut += 1
            if i in adjusted_hotspots:
                hotspot_mut += 1
    return total_mut, hotspot_mut

def adjusted_distance(pat_dist, total_mut, hotspot_mut):
    if total_mut == 0:
        return 0.0
    return pat_dist * (1 - (hotspot_mut / total_mut))

def compute_gij(aa_dict, patristic_df, hotspot_map, meta_df, verbose=False):
    gij = {}
    for i in aa_dict:
        for j in aa_dict:
            if i == j:
                continue

            # Check both cases exist in metadata and have matching serotype
            try:
                sero_i = str(meta_df.loc[i, "Serotype"]).upper()
                sero_j = str(meta_df.loc[j, "Serotype"]).upper()
                if sero_i != sero_j:
                    continue  # Skip across-serotype comparisons
            except KeyError:
                continue

            haps_i = aa_dict.get(i, [])
            haps_j = aa_dict.get(j, [])
            if not haps_i or not haps_j:
                continue

            hotspot_sites = hotspot_map.get(sero_i, [])  # use sero_i or sero_j — same here

            G_total = 0.0
            for hap_i in haps_i:
                seq_i = hap_i["seq"]
                abundance_i = hap_i["abundance"]
                header_i = hap_i["header"]
                Gi = 0.0

                for hap_j in haps_j:
                    seq_j = hap_j["seq"]
                    header_j = hap_j["header"]

                    try:
                        pat_dist = patristic_df.loc[header_i, header_j]
                    except KeyError:
                        if verbose:
                            print(f"[Missing patristic] {header_i} → {header_j}")
                        continue

                    total_mut, hotspot_mut = mutation_info(seq_i, seq_j, hotspot_sites)
                    adj_dist = adjusted_distance(pat_dist, total_mut, hotspot_mut)

                    if verbose:
                        print(f"  Pair: {header_i} → {header_j}")
                        print(f"    Patristic: {pat_dist:.4f}")
                        print(f"    Total muts: {total_mut}, Hotspot muts: {hotspot_mut}")
                        print(f"    Adj dist: {adj_dist:.4f}")

                    Gi += np.exp(adj_dist)

                G_total += Gi * abundance_i

            gij[(i, j)] = float(G_total)
    return gij

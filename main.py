from load_data import load_metadata, attach_latlon, load_aa_haplotypes
from hotspots import load_hotspots
from distances import build_distance_matrix, load_patristic_df
from compute_kij import calculate_kij
from compute_gij import compute_gij
from export import save_pij_json, export_graph
from parameters import sample_params
from weather import get_temp_over_eip
from model import compute_transmission_probs  

from pathlib import Path
from datetime import datetime
import pandas as pd
import json
def run_single_model(temp_lookup=None, run_id=None):
    print(f"[{datetime.now()}] Step 1: Load metadata")
    meta_df = load_metadata()
    dupes = meta_df.index[meta_df.index.duplicated(keep=False)]
    if not dupes.empty:
        print(f"[WARNING] Found {len(dupes)} duplicate UUID entries:")
        print(meta_df.loc[dupes].sort_index())
    assert meta_df.index.is_unique, "Duplicate UUIDs found!"

    print(f"[{datetime.now()}] Step 1a: Attach latitude/longitude")
    meta_df = attach_latlon(meta_df)
    meta_df = meta_df.dropna(subset=["Onset_of_Fever", "Blood_Collection", "lat", "lon"])
    print(f"[{datetime.now()}] Step 1b: Retained {len(meta_df)} rows with complete info.")

    print(f"[{datetime.now()}] Step 2: Load haplotypes")
    aa_dict = load_aa_haplotypes(meta_df)
    total_haps = sum(len(v) for v in aa_dict.values())
    print(f"[Step 2] Loaded haplotypes for {len(aa_dict)} UUIDs ({total_haps} total haplotypes)")

    print(f"[{datetime.now()}] Step 3: Load hotspots")
    hotspot_map = load_hotspots()
    print(f"[Step 3] Loaded hotspot positions for {len(hotspot_map)} serotypes")

    print(f"[{datetime.now()}] Step 4: Build distance matrix")
    distance_matrix = build_distance_matrix(meta_df)
    print(f"[Step 4] Distance matrix built for {len(distance_matrix)} UUIDs")

    print(f"[{datetime.now()}] Step 5: Load patristic distances")
    patristic_df = load_patristic_df()
    print(f"[Step 5] Patristic matrix shape: {patristic_df.shape}")

    print(f"[{datetime.now()}] Step 6: Sample model parameters")
    params = sample_params()
    print(f"[Step 6] Parameters: {params}")

    # Step 7: Temperature lookup (cached or fresh)
    if temp_lookup is None:
        print(f"[{datetime.now()}] Step 7: Fetch temperature from API")
        temp_lookup = {}
        success, fallback = 0, 0
        for i, row in meta_df.iterrows():
            uid = row.get("UUID", i)
            try:
                temp_lookup[uid] = get_temp_over_eip(
                    uid,
                    pd.to_datetime(row["Onset_of_Fever"]),
                    tau=0,
                    lat=row["lat"],
                    lon=row["lon"]
                )
                success += 1
            except Exception:
                temp_lookup[uid] = 28.0
                fallback += 1
        print(f"[Step 7] Temperatures fetched: {success}, defaulted: {fallback}")
    else:
        print(f"[{datetime.now()}] Step 7: Using cached temperature lookup")

    print(f"[{datetime.now()}] Step 8: Calculate Kij using mosquito dynamics")
    kij_matrix = calculate_kij(meta_df, temp_lookup, params)
    print(f"[Step 8] Kij matrix computed with {len(kij_matrix)} entries")   
    
    print(f"[{datetime.now()}] Step 9: Compute genetic distance matrix (Gij)")
    gij_matrix = compute_gij(aa_dict, patristic_df, hotspot_map, meta_df)
    print(f"[Step 9] Gij matrix computed with {len(gij_matrix)} entries")

    print(f"[{datetime.now()}] Step 10: Compute transmission matrix (Pij)")
    gamma = 1.0
    P_matrix = compute_transmission_probs(kij_matrix, gij_matrix, gamma=gamma)
    print(f"[Step 10] Pij matrix computed with {len(P_matrix)} rows")

    # Save parameters and temperature per run
    if run_id is not None:
        run_dir = Path(f"output/run_{run_id:03d}")
        run_dir.mkdir(parents=True, exist_ok=True)

        with open(run_dir / "parameters.json", "w") as f:
            json.dump(params, f, indent=2)

        with open(run_dir / "temp_lookup.json", "w") as f:
            json.dump(temp_lookup, f, indent=2)

        meta_with_temp = meta_df.copy()
        meta_with_temp["temp"] = meta_df.index.map(lambda i: temp_lookup.get(i, 28.0))
        meta_with_temp.to_csv(run_dir / "metadata_with_temp.csv")

    return P_matrix, kij_matrix, gij_matrix




def main():
    kij_matrix, gij_matrix, P_matrix = run_single_model()
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    print(f"[{datetime.now()}] Step 11: Save results")
    save_pij_json(P_matrix, str(output_dir / "P_matrix.json"))
    export_graph(P_matrix, str(output_dir / "P_matrix.graphml"))

    with open(str(output_dir / "kij_matrix.json"), "w") as f:
        json.dump(kij_matrix, f)

    with open(output_dir / "gij_matrix.json", "w") as f:
        json.dump(gij_matrix, f)
    print(f"[{datetime.now()}] Results saved to {output_dir}")

if __name__ == "__main__":
    main()

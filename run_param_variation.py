from main import run_single_model 
from weather import get_temp_over_eip
from load_data import load_metadata, attach_latlon
import pandas as pd
import numpy as np
import json
from pathlib import Path
from collections import defaultdict
from datetime import datetime

N = 100
all_pij = defaultdict(list)

# Step 1: Cache temperature once
print(f"[{datetime.now()}] Caching temperature values...")
meta_df = attach_latlon(load_metadata())
meta_df = meta_df.dropna(subset=["Onset_of_Fever", "Blood_Collection", "lat", "lon"])

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
    except Exception as e:
        temp_lookup[uid] = 28.0
        fallback += 1
print(f"[Temp caching complete] Fetched={success}, Defaulted={fallback}")

# Save temp lookup
Path("output/temperature").mkdir(parents=True, exist_ok=True)
with open("output/temperature/temp_lookup.json", "w") as f:
    json.dump(temp_lookup, f, indent=2)

# Step 2: Run model N times with saved params and outputs
for run in range(N):
    print(f"\n[Run {run + 1}/{N}] Starting model run...")
    P_matrix, _, _ = run_single_model(temp_lookup=temp_lookup, run_id=run)

    # Accumulate probabilities
    for i, row in P_matrix.items():
        for j, prob in row.items():
            all_pij[(i, j)].append(prob)

# Step 3: Summarize variation
summary = {
    f"{i}->{j}": {
        "mean": np.mean(vals),
        "std": np.std(vals),
        "q05": np.quantile(vals, 0.05),
        "q95": np.quantile(vals, 0.95)
    }
    for (i, j), vals in all_pij.items()
}

# Save summary
Path("output").mkdir(exist_ok=True)
with open("output/pij_variation_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n[Complete] Summary saved to output/pij_variation_summary.json")

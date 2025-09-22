import pandas as pd
from pathlib import Path

def load_hotspots(hotspot_dir="Hotspots"):
    """
    Load hotspot codon positions for each DENV serotype from CSV files.

    Args:
        hotspot_dir (str): Path to the directory containing DENV1.csv, DENV2.csv, etc.

    Returns:
        dict: A dictionary mapping serotypes to lists of codon positions.
    """
    hotspot_files = {
        "DENV1": Path(hotspot_dir) / "DENV1.csv",
        "DENV2": Path(hotspot_dir) / "DENV2.csv",
        "DENV3": Path(hotspot_dir) / "DENV3.csv",
    }

    hotspot_map = {}
    for sero, file_path in hotspot_files.items():
        try:
            df = pd.read_csv(file_path)
            sites = df.iloc[:, 0].dropna().astype(int).tolist()
            hotspot_map[sero.upper()] = sites
        except Exception as e:
            print(f"[Error] Failed to load {sero} hotspot file: {e}")
            hotspot_map[sero.upper()] = []

    return hotspot_map

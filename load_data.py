import pandas as pd
from Bio import SeqIO
from pathlib import Path
import os
from geopy.geocoders import Nominatim
from geopy.distance import geodesic
from time import sleep
from collections import defaultdict

manual_zip_coords = {
    100.0: (6.935584, 79.844091),
    200.0: (6.92694444, 79.84861111),
    300.0: (6.900556, 79.85333333),
    400.0: (6.888889,79.85666667),
    500.0: (6.89222222, 79.87694444),
    600.0: (6.874722,79.86055556),
    700.0: (6.90666667, 79.86333333),
    800.0: (6.91472222, 79.87777778),
    900.0: (6.93000000, 79.87777778),
    1000.0: (6.928333, 79.86416667),
    1100.0: (6.93666667, 79.84972222),
    1200.0: (6.939612, 79.862192),
    1300.0: (6.94277778, 79.85861111),
    1400.0: (6.94750000, 79.87472222),  # Updated - was significantly different
    1500.0: (6.95944444, 79.87527778),
    10107.0: (6.88864210, 79.88732270),  # Updated - minor difference
    10115.0: (6.90000000, 79.95000000),
    10118.0: (6.88266486, 79.96747806),  # Updated - minor difference
    10120.0: (6.89640000, 79.91810000),
    10150.0: (6.88970000, 79.93590000),
    10230.0: (6.85000000, 79.95000000),
    10280.0: (6.84610000, 79.92810000),
    10290.0: (6.84014470, 79.90206242),  # New
    10306.0: (6.82200000, 79.90000000),
    10350.0: (6.85324749, 79.87686789),  # Updated - significant difference
    10400.0: (6.78652835, 79.88774724),  # Updated - minor difference
    10600.0: (6.92986822, 79.89611108),  # Updated - significant difference
    10640.0: (6.93378232, 79.97978634),  # Updated - significant difference
    11300.0: (6.93000000, 79.86000000),
    11380.0: (6.94000000, 79.87000000),
    11680.0: (6.85000000, 79.91000000),
    71700.0: (7.80000000, 80.50000000),
    81000.0: (6.03000000, 80.22000000),
}

def load_metadata(path="Metadata/Metadata.csv"):
    column_names = ["UUID", "Serotype", "ZIP", "Blood_Collection", "Onset_of_Fever", "Previous_DENV"]
    df = pd.read_csv(path, names=column_names, skiprows=1)
    df["Onset_of_Fever"] = pd.to_numeric(df["Onset_of_Fever"], errors="coerce")
    df["Blood_Collection"] = pd.to_numeric(df["Blood_Collection"], errors="coerce")
    df["Onset_of_Fever"] = pd.to_datetime(df["Onset_of_Fever"], origin="1899-12-30", unit="D", errors="coerce")
    df["Blood_Collection"] = pd.to_datetime(df["Blood_Collection"], origin="1899-12-30", unit="D", errors="coerce")
    df = df.set_index("UUID")
    df.index = df.index.astype(str)
    return df


def get_latlon_for_zip(ZIP, country="Sri Lanka", cache_file="zip_cache.csv"):
    if os.path.exists(cache_file):
        cache_df = pd.read_csv(cache_file)
        if not {"ZIP", "lat", "lon"}.issubset(cache_df.columns):
            cache_df.columns = ["ZIP", "lat", "lon"]
    else:
        cache_df = pd.DataFrame(columns=["ZIP", "lat", "lon"])

    if str(ZIP) in cache_df["ZIP"].astype(str).values:
        row = cache_df[cache_df["ZIP"].astype(str) == str(ZIP)].iloc[0]
        return row["lat"], row["lon"], True

    if float(ZIP) in manual_zip_coords:
        lat, lon = manual_zip_coords[float(ZIP)]
        cache_df = pd.concat([cache_df, pd.DataFrame([{"ZIP": ZIP, "lat": lat, "lon": lon}])])
        cache_df.to_csv(cache_file, index=False)
        return lat, lon, True

    geolocator = Nominatim(user_agent="denv-model")
    try:
        location = geolocator.geocode(f"{ZIP}, {country}")
        if location:
            lat, lon = location.latitude, location.longitude
            cache_df = pd.concat([cache_df, pd.DataFrame([{"ZIP": ZIP, "lat": lat, "lon": lon}])])
            cache_df.to_csv(cache_file, index=False)
            sleep(1)
            return lat, lon, False
        else:
            return None, None, False
    except Exception:
        return None, None, False

def attach_latlon(meta_df):
    meta_df["lat"] = None
    meta_df["lon"] = None
    for idx, row in meta_df.iterrows():
        try:
            lat, lon, _ = get_latlon_for_zip(row["ZIP"])
            meta_df.at[idx, "lat"] = lat
            meta_df.at[idx, "lon"] = lon
        except:
            continue
    return meta_df

def load_aa_haplotypes(meta_df, haplo_dir="Haplotypes/"):
    hap_dict = defaultdict(list)
    uuid_set = set(map(str, meta_df.index))
    fasta_files = [f for f in os.listdir(haplo_dir) if f.endswith(".fasta")]
    for fasta_name in fasta_files:
        fasta_path = os.path.join(haplo_dir, fasta_name)
        for record in SeqIO.parse(fasta_path, "fasta"):
            header = record.id
            parts = header.split("_")
            if "consHap" not in header:
                continue
            uuid = parts[0]
            try:
                abundance = float(parts[3])
                hap_num = int(parts[2].replace("consHap", ""))
            except:
                continue
            if uuid in uuid_set:
                hap_dict[uuid].append({
                    "seq": str(record.seq),
                    "abundance": abundance,
                    "hap_id": hap_num,
                    "header": header
                })
    return dict(hap_dict)

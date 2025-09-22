import pandas as pd
from collections import defaultdict
from itertools import combinations
from geopy.distance import geodesic

def build_distance_matrix(meta_df):
    dist_matrix = defaultdict(dict)
    for idx_i, idx_j in combinations(range(len(meta_df)), 2):
        row_i = meta_df.iloc[idx_i]
        row_j = meta_df.iloc[idx_j]

        lat_i, lon_i = row_i["lat"], row_i["lon"]
        lat_j, lon_j = row_j["lat"], row_j["lon"]

        if pd.isna(lat_i) or pd.isna(lon_i) or pd.isna(lat_j) or pd.isna(lon_j):
            continue

        dist = geodesic((lat_i, lon_i), (lat_j, lon_j)).kilometers
        uuid_i = row_i.get("UUID", meta_df.index[idx_i])
        uuid_j = row_j.get("UUID", meta_df.index[idx_j])

        dist_matrix[uuid_i][uuid_j] = dist
        dist_matrix[uuid_j][uuid_i] = dist
    return dist_matrix


def load_patristic_df(path="Matrices/patristic_distances.csv"):
    return pd.read_csv(path, index_col=0)

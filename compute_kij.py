import numpy as np
import pandas as pd
from geopy.distance import geodesic
from parameters import eip_temperature, mosquito_survival, viremia_curve

D = 10  # infectious period (days)

def spatial_kernel(coord_i, coord_j, delta):
    d = geodesic(coord_i, coord_j).km
    return np.exp(-delta * d)

def calculate_kij(metadata, temp_lookup, params):
    kij = {}

    beta_h = params.get("beta_h", 0.75)
    p_mh = params.get("P_mh", 0.75)
    delta = params.get("delta", 0.001)

    for i, row_i in metadata.iterrows():
        onset_i = pd.to_datetime(row_i["Onset_of_Fever"], errors="coerce")
        if pd.isna(onset_i):
            continue
        sero_i = str(row_i["Serotype"]).upper()

        coord_i = (row_i["lat"], row_i["lon"])
        uid_i = row_i.get("UUID", i)
        T_i = temp_lookup.get(uid_i, 28.0)
        S_mosq = mosquito_survival(T_i)
        EIP_days = int(np.ceil(eip_temperature(T_i)))

        for j, row_j in metadata.iterrows():
            if i == j:
                continue

            onset_j = pd.to_datetime(row_j["Onset_of_Fever"], errors="coerce")
            if pd.isna(onset_j):
                continue

            sero_j = str(row_j["Serotype"]).upper()
            if sero_i != sero_j:
                continue  # Different serotypes → skip

            # Enforce temporal causality: j must occur after i + EIP
            min_valid_date = onset_i + pd.Timedelta(days=EIP_days)
            if onset_j < min_valid_date:
                continue

            coord_j = (row_j["lat"], row_j["lon"])
            uid_j = row_j.get("UUID", j)
            f_spatial = spatial_kernel(coord_i, coord_j, delta)

            K = sum(
                beta_h * viremia_curve(tau) * S_mosq * p_mh * f_spatial #this line
                #beta_h * viremia_curve(tau) * S_mosq * p_mh * 1.0 * f_spatial #this line
                
                for tau in range(D + 1)
            )

            kij[(uid_i, uid_j)] = K

    return kij

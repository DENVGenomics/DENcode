import pandas as pd
import requests
from datetime import timedelta

def fetch_temperature(lat, lon, date, window=10, fallback=28.0):
    try:
        start_date = (date).strftime("%Y-%m-%d")
        end_date = (date + timedelta(days=window)).strftime("%Y-%m-%d")
        url = f"https://archive-api.open-meteo.com/v1/archive?latitude={lat}&longitude={lon}&start_date={start_date}&end_date={end_date}&daily=temperature_2m_mean&timezone=UTC"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        temps = data["daily"]["temperature_2m_mean"]
        return sum(temps) / len(temps)
    except:
        return fallback

def build_temp_lookup(metadata):
    lookup = {}
    for _, row in metadata.iterrows():
        uuid = row["UUID"]
        date = pd.to_datetime(row["Onset_of_Fever"], errors="coerce")
        if pd.isna(date):
            continue
        lat = row["lat"]
        lon = row["lon"]
        temp = fetch_temperature(lat, lon, date)
        lookup[uuid] = temp
    return lookup
def get_temp_over_eip(uuid, onset_date, tau, lat, lon, eip_window=10, default_temp=28.0):
    """Fetch average temperature over the EIP window starting from onset + tau"""
    if pd.isna(onset_date):
        return default_temp
    query_date = onset_date + pd.Timedelta(days=tau)
    return fetch_temperature(lat, lon, query_date, window=eip_window, fallback=default_temp)

import pandas as pd
CSV = "/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv"
db = pd.read_csv(CSV, sep=r'\s+', usecols=["startTime", "nSamples"]).drop_duplicates()
db = db[db["nSamples"] > 0].reset_index(drop=True)
db["startTime_dt"] = pd.to_datetime(db["startTime"], errors="coerce", utc=True)
db = db.dropna(subset=["startTime_dt"]).reset_index(drop=True)
print("File 100000:", db["startTime_dt"].iloc[100000])
print("File 110079:", db["startTime_dt"].iloc[110079])

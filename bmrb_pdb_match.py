import pandas as pd
import json


bmrb_pdb_df = pd.read_csv("/Users/gaoyiting/Bioinformatik/pp1/triZOD/bmrb/bmrb_pdb.csv", sep="\t")


trizod_data = []

with open("/Users/gaoyiting/Bioinformatik/pp1/triZOD/moderate.json", "r") as f:
    for line in f:
        trizod_data.append(json.loads(line))


structure_records = []

for s in trizod_data:
    structure_records.append({
        "trizod_id": s["ID"],            # 原始 ID 用作唯一标识
        "bmrb_id": s["entryID"],    # 与 CSV 匹配用
        "sequence": s["seq"],
        "gscore": s["gscores"]           # 是个 list
    })

structure_df = pd.DataFrame(structure_records)

structure_df["bmrb_id"] = structure_df["bmrb_id"].astype(str)
bmrb_pdb_df["bmrb_id"] = bmrb_pdb_df["bmrb_id"].astype(str)

merged_df = pd.merge(structure_df, bmrb_pdb_df, on="bmrb_id")


final_df = merged_df[["pdb_id", "trizod_id", "bmrb_id", "sequence", "gscore"]]

# === 第六步：保存为新 CSV 文件 ===
final_df.to_csv("/Users/gaoyiting/Bioinformatik/pp1/triZOD/bmrb/merged_trizod_bmrb_output.csv", sep=",")

print("✅ 生成成功：merged_trizod_output.csv")

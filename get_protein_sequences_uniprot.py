import pandas as pd
import requests
import time

# Load IDs
df = pd.read_csv("dataset/final_final_filtered_dataset.csv")
uniprot_ids = df["uniprot_id"].dropna().astype(str).str.strip().unique().tolist()
print(f"Loaded {len(uniprot_ids)} UniProt IDs")

# Submit mapping job
print("Submitting ID mapping job...")
submit_response = requests.post("https://rest.uniprot.org/idmapping/run", data={
    "from": "UniProtKB_AC-ID",
    "to": "UniProtKB",
    "ids": ",".join(uniprot_ids)
})
job_id = submit_response.json()["jobId"]

# Poll for status
while True:
    status = requests.get(f"https://rest.uniprot.org/idmapping/status/{job_id}").json()
    if status.get("jobStatus") == "RUNNING":
        print("Still processing... waiting 3s")
        time.sleep(3)
    elif status.get("jobStatus") == "FINISHED" or status.get("results"):
        print("Mapping job finished.")
        break
    else:
        raise Exception("Unexpected status response")

# Fetch results
result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}?format=json"
results = requests.get(result_url).json()
mapped_entries = results.get("results", [])

# Build mapping and collect current IDs
mapped_ids = [e["to"] for e in mapped_entries]
original_to_mapped = {e["from"]: e["to"] for e in mapped_entries}

# Save mapping info
pd.DataFrame(original_to_mapped.items(), columns=["original_id", "mapped_id"]).to_csv("id_mapping_table.csv", index=False)

# Save missing (unmapped) IDs
mapped_set = set(original_to_mapped.keys())
unmapped_ids = [uid for uid in uniprot_ids if uid not in mapped_set]
with open("data/unmapped_ids.txt", "w") as f:
    for uid in unmapped_ids:
        f.write(uid + "\n")
print(f"{len(unmapped_ids)} IDs could not be mapped. Saved to 'unmapped_ids.txt'.")

# Batch-fetch FASTA sequences
def batch_iterable(iterable, size):
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]

with open("data/sequences.fasta", "w") as out_f:
    for i, batch in enumerate(batch_iterable(mapped_ids, 300)):
        print(f"Fetching batch {i+1} of {len(batch)} sequences...")
        query = " OR ".join(batch)
        url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=accession:({query})"
        r = requests.get(url)
        if r.status_code == 200:
            out_f.write(r.text)
        else:
            print(f"⚠️ Batch {i+1} failed: HTTP {r.status_code}")
        time.sleep(1)

print("✅ All available sequences saved to sequences.fasta")

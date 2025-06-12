from datasets import load_dataset

# Load the SoftDis dataset, 'clusters-all' configuration
dataset = load_dataset("CQSB/SoftDis", "clusters-all")

# Access the 'train' split
train_data = dataset["train"]
train_data.save_to_disk("rawdata/SoftDis/clusters_arrow")
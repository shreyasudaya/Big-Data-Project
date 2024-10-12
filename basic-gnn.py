import os
import torch
import torch.nn.functional as F
from torch_geometric.data import Data, DataLoader
from torch_geometric.nn import GCNConv
from Bio import SeqIO
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split

# Function to load sequences and convert them into graph data
def load_sequences(base_dir):
    data_list = []
    labels = []

    for variant_folder in os.listdir(base_dir):
        variant_path = os.path.join(base_dir, variant_folder)
        if os.path.isdir(variant_path):
            for fasta_file in os.listdir(variant_path):
                if fasta_file.endswith(".fasta"):
                    fasta_path = os.path.join(variant_path, fasta_file)
                    for record in SeqIO.parse(fasta_path, "fasta"):
                        sequence = str(record.seq).upper()
                        data = sequence_to_graph(sequence)
                        data_list.append(data)
                        labels.append(variant_folder)  # Use folder name as label

    return data_list, labels

# Function to convert a nucleotide sequence into a graph representation
def sequence_to_graph(sequence):
    # Convert bases to ASCII values as node features
    node_features = [[ord(base)] for base in sequence]
    edge_index = []
    
    for i in range(len(sequence) - 1):
        edge_index.append((i, i + 1))  # Create edges between consecutive bases

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    x = torch.tensor(node_features, dtype=torch.float)
    print("Graph Created")
    return Data(x=x, edge_index=edge_index)

# Define the GNN model
class GNNModel(torch.nn.Module):
    def __init__(self, num_features, num_classes):
        super(GNNModel, self).__init__()
        self.conv1 = GCNConv(num_features, 16)
        self.conv2 = GCNConv(16, num_classes)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        return F.log_softmax(x, dim=1)

def main():
    base_dir = "Dataset\sequences"  # Path to your dataset
    data_list, labels = load_sequences(base_dir)

    # Encode labels
    le = LabelEncoder()
    encoded_labels = le.fit_transform(labels)

    # Split data into train and test sets
    train_data, test_data, train_labels, test_labels = train_test_split(data_list, encoded_labels, test_size=0.2)

    # Create DataLoader
    train_loader = DataLoader(train_data, batch_size=32, shuffle=True)

    # Model instantiation
    num_features = 1  # Assuming single feature (ASCII value)
    num_classes = len(le.classes_)
    model = GNNModel(num_features, num_classes)

    # Optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    # Training loop
    for epoch in range(200):
        model.train()
        for data in train_loader:
            optimizer.zero_grad()
            out = model(data)
            loss = F.nll_loss(out, train_labels)
            loss.backward()
            optimizer.step()
        print(f'Epoch {epoch + 1}: Loss = {loss.item()}')

    print("Training completed.")

if __name__ == "__main__":
    main()

import torch
import numpy
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
from torch.utils.data import DataLoader

# 1. Load and Preprocess the MNIST Dataset
transform = transforms.Compose([
    transforms.ToTensor(),  # Convert images to PyTorch tensors
    transforms.Normalize((0.5,), (0.5,))  # Normalize pixel values to [-1, 1]
])

# Download training and test datasets
train_dataset = torchvision.datasets.MNIST(root='./data', train=True, download=True, transform=transform)
test_dataset = torchvision.datasets.MNIST(root='./data', train=False, download=True, transform=transform)

# Create DataLoaders
train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)

# 2. Define a Neural Network
class MNISTClassifier(nn.Module):
    def __init__(self):
        super(MNISTClassifier, self).__init__()
        self.network = nn.Sequential(
            nn.Flatten(),  # Flatten 28x28 images into a vector
            nn.Linear(28*28, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 10)  # 10 output classes (digits 0-9)
        )

    def forward(self, x):
        return self.network(x)

class MNISTClassifier(nn.Module):
    def __init__(self):
        super(MNISTClassifier, self).__init__()
        self.network = nn.Sequential(
            nn.Conv2d(1, 16, kernel_size=3, stride=1, padding=1),  # Conv layer (1 → 16 channels)
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),                 # Downsample by 2x
            nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1), # Conv layer (16 → 32 channels)
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),                 # Downsample by 2x
            nn.Flatten(),                                         # Flatten for fully connected layers
            nn.Linear(32*7*7, 128),                               # Fully connected (7x7 from pooling)
            nn.ReLU(),
            nn.Linear(128, 10)                                    # Output layer
        )

    def forward(self, x):
        return self.network(x)


# Instantiate the model
model = MNISTClassifier()

# 3. Define Loss Function and Optimizer
criterion = nn.CrossEntropyLoss()  # Suitable for multi-class classification
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 4. Train the Model
num_epochs = 5
for epoch in range(num_epochs):
    model.train()  # Set model to training mode
    for images, labels in train_loader:
        # Forward pass
        outputs = model(images)
        loss = criterion(outputs, labels)

        # Backward pass and optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}")

# 5. Evaluate the Model
model.eval()  # Set model to evaluation mode
correct = 0
total = 0
with torch.no_grad():
    for images, labels in test_loader:
        outputs = model(images)
        _, predicted = torch.max(outputs, 1)  # Get predicted class
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

accuracy = 100 * correct / total
print(f"Test Accuracy: {accuracy:.2f}%")

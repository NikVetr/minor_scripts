import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
from torchvision import datasets, transforms
import matplotlib.pyplot as plt
from tqdm import tqdm

# Device
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Hyperparameters
IMG_SIZE = 28
BATCH_SIZE = 128
EPOCHS = 5
TIMESTEPS = 1000
LEARNING_RATE = 2e-4

# Data: MNIST dataset
transform = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize((0.5,), (0.5,))
])

train_dataset = datasets.MNIST(root="./data", train=True, transform=transform, download=True)
train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)

# Beta schedule (linear for simplicity)
betas = torch.linspace(0.0001, 0.02, TIMESTEPS).to(DEVICE)
alphas = 1.0 - betas
alpha_cumprod = torch.cumprod(alphas, dim=0)
sqrt_alpha_cumprod = torch.sqrt(alpha_cumprod)
sqrt_one_minus_alpha_cumprod = torch.sqrt(1.0 - alpha_cumprod)

# Forward diffusion (add noise to images)
def forward_diffusion(x0, t):
    noise = torch.randn_like(x0)
    sqrt_alpha_t = sqrt_alpha_cumprod[t].view(-1, 1, 1, 1)
    sqrt_one_minus_alpha_t = sqrt_one_minus_alpha_cumprod[t].view(-1, 1, 1, 1)
    xt = sqrt_alpha_t * x0 + sqrt_one_minus_alpha_t * noise
    return xt, noise

# Simple U-Net (baby version)
class UNet(nn.Module):
    def __init__(self):
        super(UNet, self).__init__()
        self.down1 = nn.Sequential(
            nn.Conv2d(65, 32, 3, padding=1),  # Change input channels from 1 to 65
            nn.ReLU(),
            nn.Conv2d(32, 64, 3, stride=2, padding=1),
            nn.ReLU()
        )
        self.down2 = nn.Sequential(
            nn.Conv2d(64, 128, 3, padding=1),
            nn.ReLU(),
            nn.Conv2d(128, 128, 3, stride=2, padding=1),
            nn.ReLU()
        )
        self.up1 = nn.Sequential(
            nn.ConvTranspose2d(128, 64, 4, stride=2, padding=1),
            nn.ReLU()
        )
        self.up2 = nn.Sequential(
            nn.ConvTranspose2d(64, 32, 4, stride=2, padding=1),
            nn.ReLU(),
            nn.Conv2d(32, 1, 3, padding=1)
        )

    def forward(self, x, t):
        # Time embedding (sinusoidal)
        t_embedding = torch.sin(t.view(-1, 1).float() * torch.arange(1, 65, device=x.device).float())
        t_embedding = t_embedding.view(-1, 64, 1, 1).expand(-1, -1, x.shape[2], x.shape[3])

        # Concatenate time embedding with input image
        x = torch.cat((x, t_embedding), dim=1)  # Shape: [B, 65, H, W]

        # Forward through U-Net
        x1 = self.down1(x)
        x2 = self.down2(x1)
        x = self.up1(x2)
        x = self.up2(x + x1)

        return x

# Model, optimizer, and loss
model = UNet().to(DEVICE)
optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)

# Training loop
model.train()
for epoch in range(EPOCHS):
    total_loss = 0
    for batch_idx, (images, _) in enumerate(tqdm(train_loader, desc=f"Epoch {epoch + 1}/{EPOCHS}")):
        images = images.to(DEVICE)
        t = torch.randint(0, TIMESTEPS, (images.shape[0],), device=DEVICE)
        noisy_images, noise = forward_diffusion(images, t)
        predicted_noise = model(noisy_images, t)
        loss = F.mse_loss(predicted_noise, noise)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

        # Print every 100 batches
        if (batch_idx + 1) % 100 == 0:
            print(f"Batch {batch_idx + 1}/{len(train_loader)}, Loss: {loss.item():.4f}")

    print(f"Epoch {epoch + 1}/{EPOCHS} completed. Average Loss: {total_loss / len(train_loader):.4f}")

# Sampling from noise¸¸ZZZ
def sample(num_samples=16):
    model.eval()
    with torch.no_grad():
        x = torch.randn((num_samples, 1, IMG_SIZE, IMG_SIZE), device=DEVICE)
        for t in reversed(range(TIMESTEPS)):
            t_tensor = torch.full((num_samples,), t, device=DEVICE)
            noise_pred = model(x, t_tensor)
            alpha_t = alphas[t]
            alpha_cumprod_t = alpha_cumprod[t]
            beta_t = betas[t]

            mean = (1 / torch.sqrt(alpha_t)) * (x - beta_t / torch.sqrt(1 - alpha_cumprod_t) * noise_pred)
            if t > 0:
                noise = torch.randn_like(x)
                x = mean + torch.sqrt(beta_t) * noise
            else:
                x = mean
    return x

# Generate and display samples
samples = sample().cpu()
plt.figure(figsize=(10, 10))
for i in range(16):
    plt.subplot(4, 4, i + 1)
    plt.imshow(samples[i].squeeze(), cmap="gray")
    plt.axis("off")
plt.show()

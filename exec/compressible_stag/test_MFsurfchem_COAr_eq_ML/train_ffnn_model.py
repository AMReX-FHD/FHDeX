import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
import matplotlib.pyplot as plt
import copy

torch.set_default_dtype(torch.float64)

## read data made from temp,theta and pressure known distributions
data_path = Path('mfsurfchem_samples.txt')

rows = []
for line in data_path.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    parts = line.split()
    if len(parts) != 5:
        continue
    rows.append([float(x) for x in parts])

data = np.asarray(rows, dtype=np.float64)
pressure = data[:, 0:1]
temp = data[:, 1:2]
theta = data[:, 2:3]
rads = data[:, 3:4]
rdes = data[:, 4:5]

X = np.concatenate([pressure, temp, theta], axis=1)
Y = np.concatenate([rads, rdes], axis=1)

print('Loaded rows:', len(data))
print('Inputs [pressure, temp, theta] shape:', X.shape)
print('Outputs [rads, rdes] shape:', Y.shape)
print('First row inputs [pressure, temp, theta:', X[0])
print('First row outputs [rads, rdes]:', Y[0])

X = torch.tensor(X, dtype=torch.float64)
Y = torch.tensor(Y, dtype=torch.float64)

X_mean, X_std = X.mean(dim=0), X.std(dim=0)
Y_mean, Y_std = Y.mean(dim=0), Y.std(dim=0)

X_std = torch.where(X_std == 0, torch.tensor(1.0), X_std)
Y_std = torch.where(Y_std == 0, torch.tensor(1.0), Y_std)

X_scaled = (X - X_mean) / X_std
Y_scaled = (Y - Y_mean) / Y_std

# Split data: 70% Train, 15% Validation, 15% Test
num_samples = X.shape[0]
indices = torch.randperm(num_samples)

train_end = int(0.7 * num_samples)
val_end = int(0.85 * num_samples)
train_idx, val_idx, test_idx = indices[:train_end], indices[train_end:val_end], indices[val_end:]

X_train, Y_train = X_scaled[train_idx], Y_scaled[train_idx]
X_val, Y_val = X_scaled[val_idx], Y_scaled[val_idx]
X_test = X[test_idx] # Keep RAW test data to verify (not standardized)

# Simple Feed-Forward Neural Network via nn.Sequential
model = nn.Sequential(
    nn.Linear(3, 256),
    nn.ReLU(),
    nn.Linear(256,256),
    nn.ReLU(),
    nn.Linear(256, 2)
)

# Training loop (Train on scaled data)

# Adam Optimizer with Mean Square Error loss function
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-9)
loss_fn = nn.MSELoss()

# Lists to store the loss values for plotting
train_losses = []
val_losses = []

epochs = 2000

# Variables to track the best model
best_val_loss = float('inf')  # Start with infinity so the first epoch is always best
best_model_state = None       # To store the best weights
best_epoch = 0

# Training
for epoch in range(epochs):
    model.train()
    optimizer.zero_grad()
    train_predictions = model(X_train)
    train_loss = loss_fn(train_predictions, Y_train)
    train_loss.backward()
    optimizer.step()

    model.eval()
    with torch.no_grad():
        val_predictions = model(X_val)
        val_loss = loss_fn(val_predictions, Y_val)
    # Log the train and val loss
    current_train_loss = train_loss.item()
    current_val_loss = val_loss.item()
    train_losses.append(current_train_loss)
    val_losses.append(current_val_loss)
    # save weights of the  model with the best val loss as final
    if current_val_loss < best_val_loss:
        best_val_loss = current_val_loss
        best_epoch = epoch + 1
        best_model_state = copy.deepcopy(model.state_dict())

    if (epoch + 1) % 50 == 0:
        print(f"Epoch [{epoch+1:4d}/{epochs}], "
              f"Train Loss: {current_train_loss:.8f}, "
              f"Val Loss: {current_val_loss:.8f} "
              f"(Best: {best_val_loss:.8f} @ Ep {best_epoch})")

# Restore the best model weights
print(f"\nTraining Complete. Restoring best model from Epoch {best_epoch} with Val Loss: {best_val_loss:.8f}")
model.load_state_dict(best_model_state)

plt.figure(figsize=(8, 5))
plt.semilogy(range(1, epochs + 1), train_losses, label='Train Loss')
plt.semilogy(range(1, epochs + 1), val_losses, label='Validation Loss')
plt.axvline(x=best_epoch, color='r', linestyle='--', label=f'Best Model (Ep {best_epoch})')

plt.xlabel('Epochs')
plt.ylabel('Mean Squared Error (Scaled) [Log Scale]')
plt.title('Training and Validation Loss')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.tight_layout()
plt.show()
plt.savefig("Training_Val_loss.png")

# Make class to save model with correct standardization of the data built into the forward mechanic
class ExportedModel(nn.Module):
    def __init__(self, core, x_mean, x_std, y_mean, y_std): # core is the model and the means and std are from the training data
        super().__init__()
        self.core = core
        # store the means and std
        self.register_buffer('x_mean', x_mean.clone().detach().to(torch.float64))
        self.register_buffer('x_std', x_std.clone().detach().to(torch.float64))
        self.register_buffer('y_mean', y_mean.clone().detach().to(torch.float64))
        self.register_buffer('y_std', y_std.clone().detach().to(torch.float64))
    # define forward which will be used in AMReX/pytorch where it will be passed in the inputs and return the outputs
    def forward(self, x):
        x_scaled = (x - self.x_mean) / self.x_std # the inputs are standardized here
        y_scaled = self.core(x_scaled) # the model is called with the scaled inputs
        y_original = (y_scaled * self.y_std) + self.y_mean # the outputs are de-standardized
        return y_original # return outputs correctly scaled outputs (scaled by the mean and std of training data!)

# call model with the training mean and std
exported = ExportedModel(model, X_mean, X_std, Y_mean, Y_std)
exported.eval() # swithc model into evaluation mode
scripted = torch.jit.script(exported) # compile python code into TorchScript that will be read by AMReX/Pytorch
scripted.save('model.pt')

# Verify if model saved correctly
with torch.no_grad():
    # Load it back from disk exactly as C++ will
    loaded_model = torch.jit.load("model.pt")
    # Pass in the raw test data tensor, it should be Float64, to match with Real in AMReX
    test_predictions = loaded_model(X_test[1])
print(f"Loaded scripted model output shape: {test_predictions.shape}")
print(f"Sample Test Predictions:\n{test_predictions}")



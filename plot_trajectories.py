import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# Load CSV files
df_intersection = pd.read_csv("trajectories_intersection.csv")
df_merging = pd.read_csv("trajectories_merging.csv")
df_overtaking = pd.read_csv("trajectories_overtaking.csv")

# Function to plot trajectories
def plot_trajectories(df, scenario_name, ax):
    vehicle_ids = df["vehicle_id"].unique()
    
    for vid in vehicle_ids:
        vehicle_data = df[df["vehicle_id"] == vid]

        # Convert to numpy arrays
        x_vals = vehicle_data["x"].to_numpy()
        y_vals = vehicle_data["y"].to_numpy()
        psi_vals = vehicle_data["psi"].to_numpy()  # Vehicle heading (in radians)

        # Plot trajectory
        ax.plot(x_vals, y_vals, marker="o", label=f"{scenario_name} - Vehicle {vid}")

        # Draw the vehicle at the FIRST trajectory point (beginning)
        if len(x_vals) > 0:
            x_start, y_start, psi_start = x_vals[0], y_vals[0], psi_vals[0]

            # Vehicle dimensions
            vehicle_length = 4.5  # meters
            vehicle_width = 2.0  # meters

            # Compute bottom-left corner, considering both length & width
            x_corner = x_start - (vehicle_length / 2) * np.cos(psi_start) + (vehicle_width / 2) * np.sin(psi_start)
            y_corner = y_start - (vehicle_length / 2) * np.sin(psi_start) - (vehicle_width / 2) * np.cos(psi_start)

            # Create and add a rotated rectangle (vehicle shape)
            rect = Rectangle((x_corner, y_corner), vehicle_length, vehicle_width,
                             angle=np.degrees(psi_start), edgecolor='red', facecolor='none', lw=2)
            
            ax.add_patch(rect)  # Add vehicle shape to plot

# Create subplots (1 row, 3 columns)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 7))

# Plot each scenario on separate subplots
plot_trajectories(df_intersection, "Intersection", ax1)
plot_trajectories(df_merging, "Merging", ax2)
plot_trajectories(df_overtaking, "Overtaking", ax3)

# Configure subplot 1 (Intersection)
ax1.set_xlabel("X Position")
ax1.set_ylabel("Y Position")
ax1.set_title("Predicted Trajectories (Intersection)")
ax1.legend()
ax1.grid()
ax1.axis("equal")  # Keep proportions realistic

# Configure subplot 2 (Merging)
ax2.set_xlabel("X Position")
ax2.set_ylabel("Y Position")
ax2.set_title("Predicted Trajectories (Merging)")
ax2.legend()
ax2.grid()
ax2.axis("equal")  # Keep proportions realistic

# Configure subplot 3 (Overtaking)
ax3.set_xlabel("X Position")
ax3.set_ylabel("Y Position")
ax3.set_title("Predicted Trajectories (Overtaking)")
ax3.legend()
ax3.grid()
ax3.axis("equal")  # Keep proportions realistic

# Show the plots
plt.tight_layout()  # Adjust layout to prevent overlapping
plt.show()

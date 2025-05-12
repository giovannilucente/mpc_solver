import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# Load CSV files
df_smooth = pd.read_csv("trajectory_smooth_curve.csv")
df_lanes_smooth = pd.read_csv("lanes_smooth_curve.csv")

df_harsh = pd.read_csv("trajectory_harsh_curve.csv")
df_lanes_harsh = pd.read_csv("lanes_harsh_curve.csv")

# Function to plot rectangles along trajectory
def plot_vehicle_rectangles(df, ax):
    vehicle_length = 4.5
    vehicle_width = 2.0

    vehicle_ids = df["vehicle_id"].unique()

    for vid in vehicle_ids:
        vehicle_data = df[df["vehicle_id"] == vid]

        for _, row in vehicle_data.iterrows():
            x = row["x"]
            y = row["y"]
            psi = row["psi"]

            # Compute bottom-left corner of rectangle
            x_corner = x - (vehicle_length / 2) * np.cos(psi) + (vehicle_width / 2) * np.sin(psi)
            y_corner = y - (vehicle_length / 2) * np.sin(psi) - (vehicle_width / 2) * np.cos(psi)

            rect = Rectangle((x_corner, y_corner), vehicle_length, vehicle_width,
                             angle=np.degrees(psi), edgecolor='blue', facecolor='none', lw=1)
            ax.add_patch(rect)

# Function to plot lane centerline
def plot_centerline(df_lanes, ax):
    centerline = df_lanes[df_lanes["lane_type"] == "center"]
    ax.plot(centerline["x"].to_numpy(), centerline["y"].to_numpy(), 'k--', label="Lane centerline", linewidth=2)

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Smooth curve scenario
plot_vehicle_rectangles(df_smooth, axes[0])
plot_centerline(df_lanes_smooth, axes[0])
axes[0].set_title("Smooth Curve Scenario")
axes[0].set_xlabel("X Position")
axes[0].set_ylabel("Y Position")
axes[0].axis("equal")
axes[0].legend()
axes[0].grid()

# Harsh curve scenario
plot_vehicle_rectangles(df_harsh, axes[1])
plot_centerline(df_lanes_harsh, axes[1])
axes[1].set_title("Harsh Curve Scenario (90° Turn)")
axes[1].set_xlabel("X Position")
axes[1].set_ylabel("Y Position")
axes[1].axis("equal")
axes[1].legend()
axes[1].grid()

plt.tight_layout()
plt.show()


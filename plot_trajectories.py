import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# Load CSV files
df_intersection = pd.read_csv("trajectory_smooth_curve.csv")
df_lanes = pd.read_csv("lanes_smooth_curve.csv")

# Function to plot trajectory for one scenario
def plot_trajectory(df, scenario_name, ax):
    vehicle_ids = df["vehicle_id"].unique()

    for vid in vehicle_ids:
        vehicle_data = df[df["vehicle_id"] == vid]

        x_vals = vehicle_data["x"].to_numpy()
        y_vals = vehicle_data["y"].to_numpy()
        psi_vals = vehicle_data["psi"].to_numpy()

        ax.plot(x_vals, y_vals, marker="o", label=f"Vehicle {vid}")

        if len(x_vals) > 0:
            x_start, y_start, psi_start = x_vals[0], y_vals[0], psi_vals[0]

            vehicle_length = 4.5
            vehicle_width = 2.0

            x_corner = x_start - (vehicle_length / 2) * np.cos(psi_start) + (vehicle_width / 2) * np.sin(psi_start)
            y_corner = y_start - (vehicle_length / 2) * np.sin(psi_start) - (vehicle_width / 2) * np.cos(psi_start)

            rect = Rectangle((x_corner, y_corner), vehicle_length, vehicle_width,
                             angle=np.degrees(psi_start), edgecolor='red', facecolor='none', lw=2)
            ax.add_patch(rect)

def plot_centerline(df_lanes, ax):
    centerline = df_lanes[df_lanes["lane_type"] == "center"]
    ax.plot(centerline["x"].to_numpy(), centerline["y"].to_numpy(), 'k--', label="Lane centerline", linewidth=2)

# Create a single plot
fig, ax = plt.subplots(figsize=(10, 7))

# Plot everything
plot_trajectory(df_intersection, "Smooth curve", ax)
plot_centerline(df_lanes, ax)

# Formatting the plot
ax.set_xlabel("X Position")
ax.set_ylabel("Y Position")
ax.set_title("Planned Trajectories (Smooth curve)")
ax.legend()
ax.grid()
ax.axis("equal")

plt.tight_layout()
plt.show()
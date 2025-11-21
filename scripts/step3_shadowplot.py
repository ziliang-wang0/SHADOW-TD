import os
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import cKDTree  
import json

def make_shadow(step_b, input_path):
    start_time = time.time()
    data = np.load(input_path)
    b_vals = data["b"]
    alpha_vals = data["alpha"]
    F_vals = data["F"]

    plot_density(b_vals, alpha_vals, F_vals, step_b, save_path=plot_path)

    elapsed = time.time() - start_time
    print(f"Step one completed successfully in: {elapsed:.2f} seconds.")

def plot_density(b_vals, alpha_vals, F_vals, step_b, save_path, max_distance=0.05): 
    print("Preparing plotting data...")

    x_vals = b_vals * np.cos(alpha_vals)
    y_vals = b_vals * np.sin(alpha_vals)

    mask = (x_vals >= -xmax) & (x_vals <= xmax) & (y_vals >= -ymax) & (y_vals <= ymax)
    x_vals = x_vals[mask]
    y_vals = y_vals[mask]
    F_vals = F_vals[mask]

    print("Generating a regular grid and finding nearest neighbors...")

    grid_spacing = step_b / 1
    grid_x, grid_y = np.mgrid[-xmax:xmax:grid_spacing, -ymax:ymax:grid_spacing]
    grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
    data_points = np.column_stack((x_vals, y_vals))

    tree = cKDTree(data_points)
    dist, idx = tree.query(grid_points, distance_upper_bound=max_distance)  
    font_size = 30
    safe_F = np.zeros(len(grid_points))  
    valid_mask = np.isfinite(dist)
    safe_F[valid_mask] = F_vals[idx[valid_mask]]
    grid_F = safe_F.reshape(grid_x.shape)
    vmin = 0
    vmax = np.percentile(grid_F, 99.9)
    # Custom colormap: Black -> Red -> Yellow -> White
    cmap = LinearSegmentedColormap.from_list("custom_cmap", [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)])

    print("Plotting...")

    plt.figure(figsize=(8, 8))
    sc = plt.pcolormesh(grid_x, grid_y, grid_F, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)

    plt.gca().set_xticks([])  
    plt.gca().set_yticks([])  
    plt.gca().set_xticklabels([])  
    plt.gca().set_yticklabels([]) 
    plt.xlim([-xmax, xmax])
    plt.ylim([-ymax, ymax])
    plt.gca().set_aspect('equal', adjustable='box')

    plt.savefig(save_path, dpi=300, pad_inches=0, bbox_inches='tight')
    plt.close()
    print(f"Image saved to: {save_path}")

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(base_dir, "..", "config", "config.json")  
    with open(config_path, 'r', encoding='utf-8') as f:
        config = json.load(f)
    step_b = config["step_b"]
    xmax = config["shadow_xmax"]
    ymax = config["shadow_ymax"]
    kappa_ff = config["kappa_ff"]
    kappa_K = config["kappa_K"]
    r_max = config["r_max"]
    r_in = config["r_in"]
    psi0_deg = config["psi0_deg"]
    theta0_deg = config["theta0_deg"]
    opt_regime = config["optical_regime"].lower()
    if opt_regime == "intermediate":
        chi = config["absorption_coefficient"]
        output_opt = f"{opt_regime}_{chi:.3f}"
    else:
        output_opt = opt_regime
    
    output_dir = os.path.join(base_dir, "..", "output")
    input_name = f"flux_rmax={r_max:.1f}_optical_{output_opt}_psi0={psi0_deg:.1f}_rin={r_in:.1f}_theta0={theta0_deg:.1f}_kappaff={kappa_ff:.3f}_kappaK={kappa_K:.3f}.npz"
    input_path = os.path.join(output_dir, input_name)
    file_name = os.path.basename(input_path).replace(".npz", "")
    plot_path = os.path.join(output_dir, file_name + ".png")
    if os.path.exists(plot_path):
       print(f"Skip this step, as the file already exists: {plot_path}")
    else:
       print("Running Step 3: Shadow Plotting")
       make_shadow(step_b, input_path)

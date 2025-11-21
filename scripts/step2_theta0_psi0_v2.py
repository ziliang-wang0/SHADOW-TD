import os
import math
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from numba import njit
from numpy import sin, cos, sqrt, sign, pi
import json

# --------------------- Flux Function --------------------------
@njit
def unit_flux(r, dr, b, alpha, cos2psi, psi0, kappa_ff, kappa_K, r_in, theta0, dphi):
    r = r + (dr / 2)
    B = 1.0 - 2.0 / r   #Default metric:1/g_{rr}=-g_{tt} = B = 1 - 2M / r
    sqrt_B = sqrt(B)
    v_r = -kappa_ff * B * sqrt(2.0 / r)
    Omega = kappa_K / sqrt(r**3)
    dOmega_dr = -1.5 * kappa_K * r**(-2.5)
    Omega_in = kappa_K / sqrt(r_in**3)
    q_minus = dOmega_dr / (4.0 * pi * psi0) * (-Omega + (r_in**2 * Omega_in) / r**2)
    sin_theta = sin(theta0)
    cos_alpha = cos(alpha)
    sqrt_arg = 1.0 - (B * b**2) / r**2
    if sqrt_arg <= 0.0:
       sqrt_arg = 0.0         

    k_rabs = sqrt(sqrt_arg)
    k_r_over_k_t = (sign(dr)) / B * k_rabs
    k_phi_over_k_t = b * cos_alpha * sin_theta
    Ut_inv_sq = B - v_r**2 / B - r**2 * cos2psi * Omega**2
    U_t = 1.0 / sqrt(Ut_inv_sq)
    g = 1.0 / (U_t * (1.0 + v_r * k_r_over_k_t + Omega * k_phi_over_k_t))
    
    if b < 5.0:
        flux = q_minus * (1.0 / k_rabs) * abs(dr) * g**4  / sqrt_B
    else:
        flux = q_minus * (r**2 * dphi / b) * g**4  / sqrt_B

    return flux
@njit
def compute_absorption(r, dr, b, dphi, chi):
    r = r + (dr / 2)
    B = 1.0 - 2.0 / r #Default metric:1/g_{rr}=-g_{tt} = B = 1 - 2M / r
    sqrt_B = sqrt(B)
    sqrt_arg = 1.0 - (B * b**2) / r**2
    if sqrt_arg <= 0.0:
       sqrt_arg = 0.0         

    k_rabs = sqrt(sqrt_arg)
    
    if b < 5.0:
        dlambda =  (1.0 / k_rabs) * abs(dr) / sqrt_B
    else:
        dlambda = (r**2 * dphi / b) / sqrt_B
    absorption_ratio = math.exp(-chi * dlambda)
    return absorption_ratio

@njit
def area_flux(r, dr, b, alpha, cos2psi, psi0, kappa_ff, kappa_K, r_in, theta0, dphi):
    B = 1.0 - 2.0 / r              #Default metric:1/g_{rr}=-g_{tt} = B = 1 - 2M / r
    v_r = -kappa_ff * B * sqrt(2.0 / r)
    Omega = kappa_K / sqrt(r**3)
    dOmega_dr = -1.5 * kappa_K * r**(-2.5)
    Omega_in = kappa_K / sqrt(r_in**3)
    Q_minus = r * dOmega_dr / (4.0 * pi) * (-Omega + (r_in**2 * Omega_in) / r**2)
    sin_theta = sin(theta0)
    cos_alpha = cos(alpha)
    sqrt_arg = 1.0 - (B * b**2) / r**2
    if sqrt_arg <= 0.0:
       sqrt_arg = 0.0         

    k_rabs = sqrt(sqrt_arg)
    k_r_over_k_t = (sign(dr)) / B * k_rabs
    k_phi_over_k_t = b * cos_alpha * sin_theta
    Ut_inv_sq = B - v_r**2 / B - r**2 * cos2psi * Omega**2
    U_t = 1.0 / sqrt(Ut_inv_sq)
    g = 1.0 / (U_t * (1.0 + v_r * k_r_over_k_t + Omega * k_phi_over_k_t))
    

    flux = Q_minus  * g**4

    return flux

# --------------------- Coordinate Transform -------------------
@njit
def transform_coords_numba(r, phi, alpha, theta0):
    sin_phi = sin(phi)
    cos_phi = cos(phi)
    sin_alpha = sin(alpha)
    cos_alpha = cos(alpha)

    x = r * sin_phi * cos_alpha
    y = r * sin_phi * sin_alpha
    z = r * cos_phi

    rotation_angle = -(pi / 2 - theta0)
    cos_rot = cos(rotation_angle)
    sin_rot = sin(rotation_angle)

    x_prime = x
    y_prime = y * cos_rot - z * sin_rot
    z_prime = y * sin_rot + z * cos_rot

    return x_prime, y_prime, z_prime

@njit
def collect_valid_points(r_arr, phi_arr, dr_arr, alpha_list, theta0, psi0):
    n = len(r_arr)
    max_records = n * len(alpha_list)
    result = np.zeros((max_records, 5))  #  [alpha, r, phi, dr, cos2psi]
    count = 0

    for i in range(n):
        r = r_arr[i]
        phi = phi_arr[i]
        dr = dr_arr[i]
        for alpha in alpha_list:
            x, y, z = transform_coords_numba(r, phi, alpha, theta0)
            if (x**2 + z**2) * math.tan(psi0)**2 >= y**2:
                denom = x**2 + y**2 + z**2
                cos2psi = (x**2 + z**2) / denom
                result[count, 0] = alpha
                result[count, 1] = r
                result[count, 2] = phi
                result[count, 3] = dr
                result[count, 4] = cos2psi
                count += 1

    return result[:count]
# --------------------- Process Single b-group -----------------
def process_b_group(b_value, group_df, r_in, dalpha, theta0, psi0, kappa_ff, kappa_K, dphi, opt_regime, chi):
    n_alpha = int(round(2 * math.pi / dalpha)) + 1
    alpha_list = np.linspace(0, 2 * math.pi, n_alpha, endpoint=True)

    filtered_df = group_df[group_df["r"] >= r_in]
    r_arr = filtered_df["r"].to_numpy()
    phi_arr = filtered_df["phi"].to_numpy()
    dr_arr = filtered_df["dr"].to_numpy()

    records_array = collect_valid_points(r_arr, phi_arr, dr_arr, alpha_list, theta0, psi0)

    records = {}
    for row in records_array:
        alpha = row[0]
        if alpha not in records:
            records[alpha] = []
        records[alpha].append((row[1], row[2], row[3], row[4]))  # r, phi, dr, cos2psi

    result_rows = []
    if opt_regime == "thin":
        for alpha, data in records.items():
            data = sorted(data, key=lambda x: x[1])  # Sorted by phi
            fluxes = [unit_flux(r, dr, b_value, alpha, cos2psi, psi0, kappa_ff, kappa_K, r_in, theta0, dphi)
                      for r, phi, dr, cos2psi in data[:-1]] 
            F = sum(fluxes)
            result_rows.append((b_value, alpha, F))
    elif opt_regime == "thick":
       for alpha, data in records.items():
           if not data:
              continue
           r, phi, dr, cos2psi = min(data, key=lambda x: x[1])     
           F = area_flux(r, dr, b_value, alpha, cos2psi, psi0, kappa_ff, kappa_K, r_in, theta0, dphi)
           result_rows.append((b_value, alpha, F))
    elif opt_regime == "intermediate":
        for alpha, data in records.items():
           data = sorted(data, key=lambda x: x[1], reverse=True)  # Sort by phi from largest to smallest
           F = 0.0  
           for r, phi, dr, cos2psi in data:
               flux = unit_flux(r, dr, b_value, alpha, cos2psi, psi0,
                     kappa_ff, kappa_K, r_in, theta0, dphi)
               ab_ratio = compute_absorption(r, dr, b_value, dphi, chi)
               F = F * ab_ratio + flux
           result_rows.append((b_value, alpha, F))
    else:
        raise ValueError(f"Unknown optical regime: '{opt_regime}'. " "Expected one of: 'thin', 'thick', or 'intermediate'.")
    return result_rows

# --------------------- Chunk Worker ---------------------------
def process_b_range(min_b, max_b, input_file, temp_file, r_in, dalpha, theta0, psi0, kappa_ff, kappa_K, dphi, step_b, r_max, opt_regime, chi):
    results = []

    raw_array = np.load(input_file)  
    df = pd.DataFrame(raw_array, columns=["b", "r", "phi", "dr"])
    chunk = df[(df["b"] >= min_b) & (df["b"] < max_b)]
    for b_val, group in chunk.groupby("b"):
        results.extend(process_b_group(b_val, group, r_in, dalpha, theta0, psi0, kappa_ff, kappa_K, dphi, opt_regime, chi))

    out_df = pd.DataFrame(results, columns=["b", "alpha", "F"])
    np.save(temp_file, out_df.to_numpy())

# --------------------- Main ------------------------
def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(base_dir, "..", "config", "config.json")  
    with open(config_path, 'r', encoding='utf-8') as f:
        config = json.load(f)

    kappa_ff = config["kappa_ff"]
    kappa_K = config["kappa_K"]
    r_in = config["r_in"]
    r_max = config["r_max"]
    psi0_deg = config["psi0_deg"]
    theta0_deg = config["theta0_deg"]
    dalpha = config["dalpha"]
    dphi = config["dphi"]
    step_b = config["step_b"]
    theta0 =  math.radians(theta0_deg) 
    psi0 = math.radians(psi0_deg)
    opt_regime = config["optical_regime"].lower()
    chi = config["absorption_coefficient"]
    if opt_regime == "intermediate":
        output_opt = f"{opt_regime}_{chi:.3f}"
    else:
        output_opt = opt_regime
    output_dir = os.path.join(base_dir, "..", "output")
    save_path_name = f"flux_rmax={r_max:.1f}_optical_{output_opt}_psi0={psi0_deg:.1f}_rin={r_in:.1f}_theta0={theta0_deg:.1f}_kappaff={kappa_ff:.3f}_kappaK={kappa_K:.3f}.npz"
    npz_path = os.path.join(output_dir, save_path_name)
    max_workers = max(1, os.cpu_count() // 4)
    if os.path.exists(npz_path):
       print(f"Skip this step, as the file already exists: {npz_path}")
    else:
       print("Running Step 2: Flux Calculation")
       start_time = time.time()
       input_name = f"geodesics_rmax={r_max:.1f}_step_b={step_b:.3f}_dphi={dphi:.3f}_geothick.npy"
       input_file = os.path.join(output_dir, input_name)
       temp_folder = os.path.join(output_dir, "temp_chunks")
       os.makedirs(temp_folder, exist_ok=True)

    # Clean up old temp files
       for f in os.listdir(temp_folder):
           os.remove(os.path.join(temp_folder, f))

       ranges = [(i * 0.1, (i + 1) * 0.1) for i in range(250)]
       tasks = []
       futures_to_range = {}
       with ProcessPoolExecutor(max_workers=max_workers) as executor:
           for i, (min_b, max_b) in enumerate(ranges):
               temp_file = os.path.join(temp_folder, f"chunk_{i:03d}.npy")
               future = executor.submit(process_b_range, min_b, max_b, input_file, temp_file, r_in, dalpha, theta0, psi0, kappa_ff, kappa_K, dphi, step_b, r_max, opt_regime, chi)
               tasks.append(future)
               futures_to_range[future] = (min_b, max_b)
     
           for future in tqdm(as_completed(tasks), total=len(tasks), desc="Processing range b"):
               try:
                   result = future.result()
               except Exception as e:
                   min_b, max_b = futures_to_range[future]
                   print(f"Processing range [{min_b}, {max_b}) failed: {e}")

       valid_arrays = []

       for f in sorted(os.listdir(temp_folder)):
           if f.endswith(".npy"):
               path = os.path.join(temp_folder, f)
               try:
                   array = np.load(path, allow_pickle=True)
                   if array.size == 0:
                    continue
                   if array.shape[1] != 3:
                       print(f"{path} Dimension mismatch, skip.")
                       continue
                   valid_arrays.append(array)
               except Exception as e:
                   print(f"Reading {path} failed: {e}")
                   continue

       if valid_arrays:
    # Merge all arrays and sort
           merged_array = np.vstack(valid_arrays)
           sorted_array = merged_array[np.lexsort((merged_array[:,1], merged_array[:,0]))]  # Sort by b, alpha

     #Split into three columns
           b_all = sorted_array[:, 0]
           alpha_all = sorted_array[:, 1]
           F_all = sorted_array[:, 2]

    # save
           npz_path = os.path.join(output_dir, save_path_name)
           np.savez_compressed(npz_path, b=b_all, alpha=alpha_all, F=F_all)

           print(f"Data saved to:{npz_path}")
       else:
           print("No valid data is available for merging.")

       print(f"Step two completed successfully in:{time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()

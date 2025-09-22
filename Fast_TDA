


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.linalg
import os
from gtda.time_series import SingleTakensEmbedding

#  Ellipse functions
# --------------------

def rotate_ellipse(A, B, C, D, E, F):
    theta = 0.5 * np.arctan2(B, A - C)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    A_prime = A * cos_t**2 + B * cos_t * sin_t + C * sin_t**2
    B_prime = 0.0
    C_prime = A * sin_t**2 - B * cos_t * sin_t + C * cos_t**2

    D_prime = D * cos_t + E * sin_t
    E_prime = -D * sin_t + E * cos_t
    F_prime = F

    return A_prime, B_prime, C_prime, D_prime, E_prime, F_prime, theta



def ellipse_to_standard_form(A, B, C, D, E, F):
    
#  B == 0 
# Ellipse Center h,k
    
    h = -D / (2*A)
    k = -E / (2*C)

    F_new = F - (D**2)/(4*A) - (E**2)/(4*C)

    if not (A > 0 and C > 0 and F_new < 0):
        raise ValueError("Imaginary ellipse")

# a should be major axis
    
    a = np.sqrt(-F_new / A)
    b = np.sqrt(-F_new / C)

    if a < b:
        a, b = b, a

    return (h, k), a, b

def max_consecutive_distance(X, Y):
    dists = np.sqrt(np.diff(X)**2 + np.diff(Y)**2)
    return np.max(dists)


# Embeddding Theorem Parameters 
# Giotto-TDA package
# ------------------------------

window_size = 1320
step_size = 10
time_delays = 30 
dimension = 2  
stride = 1  
n_jobs = -1  # CPU



#  Load Data 
# ----------

data_file = 'Temp_21_output_time.csv'
data = pd.read_csv(data_file, header=0)
dt = data.iloc[2,0] - data.iloc[1,0]
input_time = data.iloc[:,0]
input_f = data.iloc[:,1]
data_length = len(data)

#  Counters 
failed_count = 0
total_count = 0

# saving features
# --------------------------------

ellipse_features = {
    "center_x": [],
    "center_y": [],
    "semi_major": [],
    "semi_minor": [],
    "angle": [],
    "max_consec_dist": []
}
time_axis = []

#  Main Loop 
# -----------

for j in range(window_size, data_length, step_size):
    total_count += 1
    current_window = data.iloc[j-window_size:j]
    data_raw = current_window.iloc[:,1].values.flatten()

    embedding = SingleTakensEmbedding('fixed', time_delays, dimension, stride, n_jobs)
    point_cloud = embedding.fit_transform(data_raw)
    X, Y = point_cloud[:,0], point_cloud[:,1]

#  A + C = 1
#----------    
    G = np.vstack([X**2 - Y**2, X*Y, X, Y, np.ones_like(X)]).T
    t = -(Y**2)

    try:
        U, S, Vt = np.linalg.svd(G, full_matrices=False)
    except np.linalg.LinAlgError:
        failed_count += 1
        continue

# singular values
#---------------    
    rcond = np.finfo(float).eps * max(G.shape)
    S_inv = np.diag(np.where(S > rcond * (S[0] if S.size else 1.0), 1.0/S, 0.0))

    A_est, B_est, D_coef, E_est, F_est = Vt.T @ S_inv @ (U.T @ t)

    A, B = A_est, B_est
    C = 1.0 - A_est
    D, E, F = D_coef, E_est, F_est

    if 4*A*C - B**2 <= 0:
        failed_count += 1
        continue

#  Convert to standard ellipse form
# ---------------------------------
    try:
        A_rot, B_rot, C_rot, D_rot, E_rot, F_rot, theta = rotate_ellipse(A, B, C, D, E, F)
        center, semi_major_axis, semi_minor_axis = ellipse_to_standard_form(
            A_rot, B_rot, C_rot, D_rot, E_rot, F_rot
        )
    except ValueError:
        failed_count += 1
        continue

    x0, y0 = center

#  save features
#---------------
    
    ellipse_features["center_x"].append(x0)
    ellipse_features["center_y"].append(y0)
    ellipse_features["semi_major"].append(semi_major_axis)
    ellipse_features["semi_minor"].append(semi_minor_axis)
    ellipse_features["angle"].append(theta)
    ellipse_features["max_consec_dist"].append(max_consecutive_distance(X, Y))

    time_axis.append(input_time[j - 1])

# --- Pack results into a DataFrame and save ---
features_df = pd.DataFrame({
    "time": time_axis,
    "center_x": ellipse_features["center_x"],
    "center_y": ellipse_features["center_y"],
    "semi_major": ellipse_features["semi_major"],
    "semi_minor": ellipse_features["semi_minor"],
    "angle": ellipse_features["angle"],
    "max_consec_dist": ellipse_features["max_consec_dist"],
})

features_df.to_csv("ellipse_features.csv", index=False)

print(f"Total attempts: {total_count}")
print(f"Failed fits: {failed_count}")
print(f"Successful fits: {total_count - failed_count}")


# Final Feature Plot 
# ------------------

plt.rcParams["font.family"] = "Times New Roman"
plt.figure(figsize=(12, 8))
keys = ["center_x", "center_y", "semi_major", "semi_minor", "angle", "max_consec_dist"]

for idx, key in enumerate(keys):
    plt.subplot(3, 2, idx+1)
    plt.plot(time_axis, ellipse_features[key])
    plt.title(key)
    plt.axvline(x=0.020, color='red', linestyle='--') # impact time

plt.tight_layout()
plt.show()


# pick one: "center_x", "center_y", "semi_major", "semi_minor", "angle", "max_consec_dist"
# ---------------------------------------------------------------------------------------

feature_key = "semi_major"   
feature = ellipse_features[feature_key]

# Define threshold (approximately)
# -------------------------------
threshold = 100

plt.figure(figsize=(12, 5))
plt.plot(input_time, input_f, color='black', label='Time series')

for i in range(len(time_axis) - 1):
    if feature[i] >= threshold:
        plt.axvspan(time_axis[i], time_axis[i+1], facecolor='lightblue', alpha=0.3)
    else:
        plt.axvspan(time_axis[i], time_axis[i+1], facecolor='pink', alpha=0.3)

plt.xlabel("Time (s)")
plt.ylabel("Signal")
plt.title(f"Time Series with Shading by Feature: {feature_key}")
plt.legend()
plt.show()


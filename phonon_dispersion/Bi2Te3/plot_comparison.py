import matplotlib
matplotlib.use('Agg')  # for cluster (no GUI)

import yaml
import numpy as np
import matplotlib.pyplot as plt
import os
import sys


# -------------------------------
# Debug helper
# -------------------------------
def debug_print(msg):
    print(f"[DEBUG] {msg}")


# -------------------------------
# Load YAML safely
# -------------------------------
def load_band_yaml(filename):
    debug_print(f"Checking file: {filename}")

    if not os.path.exists(filename):
        print(f"[ERROR] File not found: {filename}")
        sys.exit(1)

    debug_print(f"Loading YAML: {filename}")
    with open(filename, 'r') as f:
        try:
            data = yaml.safe_load(f)
        except Exception as e:
            print(f"[ERROR] YAML parsing failed: {e}")
            sys.exit(1)

    # Check keys
    if 'phonon' not in data:
        print(f"[ERROR] 'phonon' key missing in {filename}")
        sys.exit(1)

    phonon_data = data['phonon']
    debug_print(f"Number of q-points: {len(phonon_data)}")

    distances = []
    frequencies = []

    for i, point in enumerate(phonon_data):
        if 'distance' not in point:
            print(f"[ERROR] Missing 'distance' at q-point {i}")
            sys.exit(1)

        if 'band' not in point:
            print(f"[ERROR] Missing 'band' at q-point {i}")
            sys.exit(1)

        distances.append(point['distance'])

        freqs = []
        for j, band in enumerate(point['band']):
            if 'frequency' not in band:
                print(f"[ERROR] Missing 'frequency' at q-point {i}, band {j}")
                sys.exit(1)
            freqs.append(band['frequency'])

        frequencies.append(freqs)

    distances = np.array(distances)
    frequencies = np.array(frequencies)

    debug_print(f"{filename} loaded successfully")
    debug_print(f"Shape: {frequencies.shape}")

    return distances, frequencies


# -------------------------------
# MAIN
# -------------------------------
debug_print("Starting script...")

# Load both files
dist_dft, freq_dft = load_band_yaml('band_dft.yaml')
dist_nep, freq_nep = load_band_yaml('band_nep.yaml')

# -------------------------------
# Consistency checks
# -------------------------------
debug_print("Checking consistency...")

if len(dist_dft) != len(dist_nep):
    print("[WARNING] Different number of q-points!")
    print(f"DFT: {len(dist_dft)}, NEP: {len(dist_nep)}")

if freq_dft.shape[1] != freq_nep.shape[1]:
    print("[WARNING] Different number of bands!")
    print(f"DFT bands: {freq_dft.shape[1]}, NEP bands: {freq_nep.shape[1]}")

# Check for NaN or weird values
if np.isnan(freq_dft).any():
    print("[WARNING] NaN values found in DFT frequencies")

if np.isnan(freq_nep).any():
    print("[WARNING] NaN values found in NEP frequencies")

# Check negative frequencies
if (freq_dft < 0).any():
    print("[INFO] DFT has imaginary modes (negative frequencies)")

if (freq_nep < 0).any():
    print("[INFO] NEP has imaginary modes (negative frequencies)")


# -------------------------------
# Plot
# -------------------------------
debug_print("Plotting...")

plt.figure(figsize=(6, 4))

# DFT
for i in range(freq_dft.shape[1]):
    plt.plot(dist_dft, freq_dft[:, i], 'k-', lw=1.5,
             label='DFT' if i == 0 else "")

# NEP
for i in range(freq_nep.shape[1]):
    plt.plot(dist_nep, freq_nep[:, i], 'r--', lw=1.2,
             label='NEP' if i == 0 else "")

plt.xlabel('Wave Vector')
plt.ylabel('Frequency (THz)')
plt.title('Phonon Band: DFT vs NEP')
plt.legend()

# -------------------------------
# Save
# -------------------------------
output_file = "band_comparison.png"
plt.tight_layout()
plt.savefig(output_file, dpi=300)

debug_print(f"Plot saved as {output_file}")


# -------------------------------
# Optional: RMS error calculation
# -------------------------------
debug_print("Computing RMS error...")

try:
    min_len = min(len(dist_dft), len(dist_nep))
    error = freq_dft[:min_len] - freq_nep[:min_len]
    rms = np.sqrt(np.mean(error**2))
    print(f"[RESULT] RMS Error (THz): {rms:.4f}")
except Exception as e:
    print(f"[WARNING] RMS calculation failed: {e}")


debug_print("Script completed successfully 🚀")

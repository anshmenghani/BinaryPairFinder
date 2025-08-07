import pandas as pd
import tkinter as tk
import time
from tkinter import filedialog, messagebox, scrolledtext
import networkx as nx
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u


df = pd.DataFrame()
clusters = []
wds_coords = None
log_messages = []

def log(msg):
    log_messages.append(msg)
    results_text.insert(tk.END, msg + "\n")
    results_text.see(tk.END)

def load_csv():
    global df, clusters
    file_path = filedialog.askopenfilename(
        title="Select Gaia CSV File",
        filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
    )
    if file_path:
        try:
            df = pd.read_csv(file_path, dtype={'source_id': str})
            if 'parallax' not in df.columns:
                messagebox.showerror("Error", "CSV file must contain 'parallax' column.")
                return
            df['distance_pc'] = 1000.0 / df['parallax']
            clusters = []
            log(f"\nLoaded file: {file_path}")
            log(f"Entries loaded: {len(df)}")
        except Exception as e:
            messagebox.showerror("Error loading CSV", str(e))

def format_entry(index):
    try:
        sid = df.iloc[index]['source_id']
        if pd.notna(sid) and sid != '':
            return f"{index} (source_id={sid})"
    except:
        pass
    return f"{index}"

def hms_to_degrees(hms_str):
    try:
        parts = list(map(float, hms_str.strip().split(':')))
        if len(parts) != 3:
            raise ValueError("Invalid HMS format. Use HH:MM:SS")
        h, m, s = parts
        return (h + m / 60 + s / 3600) * 15
    except Exception:
        raise ValueError("Failed to parse RA HMS input (HH:MM:SS)")

def dms_to_degrees(dms_str):
    try:
        parts = list(map(float, dms_str.strip().split(':')))
        if len(parts) != 3:
            raise ValueError("Invalid DMS format. Use DD:MM:SS")
        d, m, s = parts
        sign = -1 if d < 0 else 1
        return sign * (abs(d) + m / 60 + s / 3600)
    except Exception:
        raise ValueError("Failed to parse Dec DMS input (DD:MM:SS)")

def find_matches(min_ra_sep, max_ra_sep, min_dec_sep, max_dec_sep, ruwe_thresh, dist_thresh, min_ruwe, max_ruwe, timed_rows=False):
    if timed_rows:
        df_sampled = df.sample(n=timed_rows, random_state=42).reset_index(drop=True)
    else:
        df_sampled = df.copy()
    G = nx.Graph()
    for i, row_i in df_sampled.iterrows():
        # Filter stars outside RUWE range
        if row_i['ruwe'] < min_ruwe or row_i['ruwe'] > max_ruwe:
            continue
        G.add_node(i)
        for j in range(i + 1, len(df_sampled)):
            row_j = df_sampled.iloc[j]
            if row_j['ruwe'] < min_ruwe or row_j['ruwe'] > max_ruwe:
                continue
            ra_diff = abs(row_i['ra'] - row_j['ra'])
            dec_diff = abs(row_i['dec'] - row_j['dec'])
            ruwe_diff = abs(row_i['ruwe'] - row_j['ruwe'])
            dist_diff = abs(row_i['distance_pc'] - row_j['distance_pc'])
            if (min_ra_sep <= ra_diff <= max_ra_sep and
                min_dec_sep <= dec_diff <= max_dec_sep and
                ruwe_diff <= ruwe_thresh and
                dist_diff <= dist_thresh):
                G.add_edge(i, j)

    grouped = [list(component) for component in nx.connected_components(G) if len(component) >= 2]
    return grouped

def run_matching():
    global clusters
    if df.empty:
        messagebox.showwarning("No Data", "Please load a CSV file first.")
        return
    try:
        min_ra_sep = hms_to_degrees(min_ra_hms_entry.get())
        max_ra_sep = hms_to_degrees(max_ra_hms_entry.get())
        min_dec_sep = dms_to_degrees(min_dec_dms_entry.get())
        max_dec_sep = dms_to_degrees(max_dec_dms_entry.get())
        ruwe_thresh = float(ruwe_entry.get())
        dist_thresh = float(dist_entry.get())
        min_ruwe = float(min_ruwe_entry.get())
        max_ruwe = float(max_ruwe_entry.get())

        # Input validation for ranges
        if min_ra_sep > max_ra_sep:
            messagebox.showerror("Input Error", "Minimum RA separation cannot be greater than maximum RA separation.")
            return
        if min_dec_sep > max_dec_sep:
            messagebox.showerror("Input Error", "Minimum Dec separation cannot be greater than maximum Dec separation.")
            return
        if min_ruwe > max_ruwe:
            messagebox.showerror("Input Error", "Minimum RUWE cannot be greater than maximum RUWE.")
            return

        log("\nRunning match algorithm...")
        clusters = find_matches(min_ra_sep, max_ra_sep, min_dec_sep, max_dec_sep, ruwe_thresh, dist_thresh, min_ruwe, max_ruwe)
        log(f"Matching complete. Found {len(clusters)} group(s).")

        if not clusters:
            log("No matches found.")
        else:
            for idx, group in enumerate(clusters):
                entries = [format_entry(i) for i in group]
                log(f"Group {idx + 1} | {len(group)} entries | index, (source_id): {entries}")
    except Exception as e:
        messagebox.showerror("Error", str(e))

def export_to_csv():
    global clusters
    if not clusters:
        messagebox.showinfo("Nothing to export", "No matching groups found yet.")
        return
    try:
        matched_rows = []
        for group_id, group in enumerate(clusters, 1):
            for idx in group:
                row = df.iloc[idx].copy()
                row['group_id'] = group_id
                matched_rows.append(row)
        result_df = pd.DataFrame(matched_rows)

        if 'source_id' in result_df.columns:
            result_df['source_id'] = result_df['source_id'].fillna('')

        save_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if save_path:
            result_df.to_csv(save_path, index=False)
            log(f"\nExported results to: {save_path}")
            messagebox.showinfo("Success", f"Exported to {save_path}")
    except Exception as e:
        messagebox.showerror("Export Error", str(e))

def load_wds_catalog():
    file_path = filedialog.askopenfilename(
        title="Select WDS Text File",
        filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
    )
    if not file_path:
        return None

    coords = []
    skipped = 0
    with open(file_path, 'r') as f:
        for line in f:
            if len(line) < 18:
                skipped += 1
                continue
            coord_str = line[-18:].strip()

            try:
                # RA: first 8 characters (HHMMSS.ss)
                ra_str = coord_str[0:8]
                # Dec: rest (±DDMMSS.s)
                dec_str = coord_str[8:]

                if '.' not in ra_str or '.' not in dec_str:
                    skipped += 1
                    continue

                # Convert RA to hours
                ra_h = float(ra_str[0:2])
                ra_m = float(ra_str[2:4])
                ra_s = float(ra_str[4:])
                ra_hours = ra_h + ra_m/60 + ra_s/3600

                # Convert Dec to degrees
                sign = -1 if dec_str[0] == '-' else 1
                dec_d = float(dec_str[1:3])
                dec_m = float(dec_str[3:5])
                dec_s = float(dec_str[5:])
                dec_deg = sign * (dec_d + dec_m/60 + dec_s/3600)

                coords.append(SkyCoord(ra=ra_hours*u.hourangle, dec=dec_deg*u.deg))
            except Exception:
                skipped += 1
                continue

    log(f"\nLoaded WDS catalog with {len(coords)} valid entries "
        f"(skipped {skipped} malformed lines) from: {file_path}")
    return SkyCoord(coords)

def compare_with_wds():
    global df, wds_coords
    if df.empty:
        messagebox.showwarning("Warning", "Please load a Gaia CSV file first.")
        return

    if wds_coords is None:
        wds_coords = load_wds_catalog()
        if wds_coords is None or len(wds_coords) == 0:
            messagebox.showerror("Error", "Failed to load WDS catalog or catalog is empty.")
            return

    try:
        gaia_coords = SkyCoord(ra=df['ra'].values * u.degree, dec=df['dec'].values * u.degree)
        idx, d2d, _ = gaia_coords.match_to_catalog_sky(wds_coords)
        matches = d2d < 1 * u.arcsec
        df['WDS_match'] = matches
        matched_count = matches.sum()

        log(f"\nWDS comparison complete. Matches found: {matched_count}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed during coordinate matching:\n{e}")

# GUI SETUP
root = tk.Tk()

try:
    icon_img = tk.PhotoImage(file="BinaryStarLogo.ico")
    root.iconphoto(True, icon_img)
except Exception as e:
    # No results_text yet, so just print
    print(f"Could not load icon: {e}")

root.title("Binary Star Finder")

title_label = tk.Label(root, text="Binary Star Finder", font=("Helvetica", 18, "bold"))
title_label.grid(row=0, column=0, columnspan=6, pady=(10, 0))
author_label = tk.Label(root, text="By: Ansh Menghani, ansh.menghani@gmail.com", font=("Helvetica", 10, "italic"))
author_label.grid(row=1, column=0, columnspan=6, pady=(0, 10))

tk.Button(root, text="Select Gaia CSV File", command=load_csv).grid(row=2, column=0, columnspan=6, pady=5)

# Min and Max RA Thresholds (HH:MM:SS)
tk.Label(root, text="Min RA Separation (HH:MM:SS):").grid(row=3, column=0)
min_ra_hms_entry = tk.Entry(root)
min_ra_hms_entry.insert(0, "00:00:00")
min_ra_hms_entry.grid(row=3, column=1)

tk.Label(root, text="Max RA Separation (HH:MM:SS):").grid(row=3, column=2)
max_ra_hms_entry = tk.Entry(root)
max_ra_hms_entry.insert(0, "00:06:00")  # approx 0.1 deg
max_ra_hms_entry.grid(row=3, column=3)

# Min and Max Dec Thresholds (DD:MM:SS)
tk.Label(root, text="Min Dec Separation (DD:MM:SS):").grid(row=4, column=0)
min_dec_dms_entry = tk.Entry(root)
min_dec_dms_entry.insert(0, "00:00:00")
min_dec_dms_entry.grid(row=4, column=1)

tk.Label(root, text="Max Dec Separation (DD:MM:SS):").grid(row=4, column=2)
max_dec_dms_entry = tk.Entry(root)
max_dec_dms_entry.insert(0, "00:06:00")  # approx 0.1 deg
max_dec_dms_entry.grid(row=4, column=3)

# RUWE inputs
tk.Label(root, text="RUWE Difference Threshold:").grid(row=5, column=0)
ruwe_entry = tk.Entry(root)
ruwe_entry.insert(0, "0.5")
ruwe_entry.grid(row=5, column=1)

tk.Label(root, text="RUWE Lower Bound:").grid(row=5, column=2)
min_ruwe_entry = tk.Entry(root)
min_ruwe_entry.insert(0, "1.2")
min_ruwe_entry.grid(row=5, column=3)

tk.Label(root, text="RUWE Upper Bound:").grid(row=6, column=0)
max_ruwe_entry = tk.Entry(root)
max_ruwe_entry.insert(0, "2.0")
max_ruwe_entry.grid(row=6, column=1)

# Distance threshold
tk.Label(root, text="Distance Threshold (pc):").grid(row=6, column=2)
dist_entry = tk.Entry(root)
dist_entry.insert(0, "2")
dist_entry.grid(row=6, column=3)

tk.Button(root, text="Find Matches", command=run_matching).grid(row=7, column=0, columnspan=6, pady=5)
tk.Button(root, text="Export to CSV", command=export_to_csv).grid(row=8, column=0, columnspan=6)
tk.Button(root, text="Compare with WDS Catalog", command=compare_with_wds).grid(row=9, column=0, columnspan=6, pady=5)

results_text = scrolledtext.ScrolledText(root, width=120, height=20)
results_text.grid(row=10, column=0, columnspan=6)

copyright_label = tk.Label(
    root,
    text="© 2025 Ansh Menghani. Binary Star Finder. All rights reserved.",
    font=("Helvetica", 8),
    fg="gray"
)
copyright_label.grid(row=11, column=0, columnspan=6, pady=(5, 10))

def show_instructions():
    instr_win = tk.Toplevel(root)
    instr_win.title("Instructions & Definitions")
    instr_win.geometry("700x500")

    text_area = tk.Text(instr_win, wrap=tk.WORD, padx=10, pady=10)
    scrollbar = tk.Scrollbar(instr_win, command=text_area.yview)
    text_area.configure(yscrollcommand=scrollbar.set)

    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    text_area.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    instructions = """
Binary Star Finder Instructions & Definitions

1. Parameter Definitions:
-------------------------
- Min RA Separation (HH:MM:SS):
  Minimum allowed difference in Right Ascension (in degrees) between stars to consider them as a pair.

- Max RA Separation (HH:MM:SS):
  Maximum allowed difference in Right Ascension (in degrees) between stars to consider them as a pair.

- Min Dec Separation (DD:MM:SS):
  Minimum allowed difference in Declination (in degrees) between stars for pairing.

- Max Dec Separation (DD:MM:SS):
  Maximum allowed difference in Declination (in degrees) between stars for pairing.

- RUWE Difference Threshold:
  Maximum allowed difference in RUWE (Renormalized Unit Weight Error) values. RUWE indicates data quality; similar RUWE values suggest similar data reliability. RUWE also indicated how well the astrometric solution fits the data. A higher RUWE (>~1.2) indicates a poorer fit to the single-star solution, potentially indicating binary stars or other light sources.

- Minimum RUWE:
  Minimum RUWE value to filter stars with poor astrometric fits. Only stars with RUWE >= this value are considered.

- Maximum RUWE:
  Maximum RUWE value to exclude stars with poor astrometric fits above this limit.

- Distance Threshold (pc):
  Maximum difference in distance (parsecs) between stars. Ensures physical proximity in space.

2. Reasoning for Choosing Parameters:
------------------------------------
These thresholds help identify candidate binary stars by filtering for stars that are close together on the sky, at similar distances, and with similar data quality (RUWE). Tight thresholds reduce false positives but might miss wider binaries.

3. Gaia Catalog:
----------------
Gaia is a space observatory that provides precise astrometric data (positions, parallaxes, proper motions) for over a billion stars. Your input CSV should contain at least:

- source_id (string): Unique star identifier.
- ra (float): Right Ascension in degrees.
- dec (float): Declination in degrees.
- parallax (float): Parallax in milliarcseconds (mas), used to compute distance.
- ruwe (float): Astrometric goodness-of-fit indicator.

Example CSV header (minimum required columns):
source_id,ra,dec,parallax,ruwe
1234567890123456789,150.1234,-35.5678,10.5,1.25

4. WDS Catalog:
---------------
Washington Double Star Catalog (WDS) lists known binary stars with precise coordinates. The app compares Gaia stars to WDS entries to identify known binaries.

Your WDS catalog should be a text file with star coordinates at the end of each line, for example:

... 12 34 56.7 +12 34 56.7

Coordinates are in the format HH MM SS.S (RA) and ±DD MM SS.S (Dec). If you include headers, they will be counted as malformed lines and skipped during WDS-Gaia matching.

5. File Formats:
----------------
- Gaia CSV: Comma-separated values with headers, UTF-8 encoded, containing required columns.
- WDS TXT: Plain text file with fixed-width lines ending with coordinate strings parseable by astropy's SkyCoord.

---

Feel free to adjust thresholds to balance sensitivity and precision. Contact Ansh Menghani at ansh.menghani@gmail.com for support.
"""

    text_area.insert(tk.END, instructions)
    text_area.configure(state=tk.DISABLED)

instr_button = tk.Button(root, text="Instructions & Definitions", command=show_instructions)
instr_button.grid(row=12, column=0, columnspan=6, pady=(5, 10))

root.mainloop()

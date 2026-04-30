"""
Visualization: Overlay Visium HD Capture Area on H&E Tissue Image

Loads the sketching results CSV and the Visium spatial data
to overlay the sequenced capture area on the full tissue scan.

Also shows a side-by-side of full tissue vs. sketch subset.

For Jupyter: copy the NOTEBOOK block below into a cell.
For CLI:     python visualize_sketch.py
"""

# ══════════════════════════════════════════════════════════════════════
# NOTEBOOK VERSION — copy everything below into a Jupyter cell
# ══════════════════════════════════════════════════════════════════════
#
# import pandas as pd
# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# import json
#
# # ── Paths (adjust resolution folder if needed) ──
# CSV_PATH          = "../bin/visium_real_data_results.csv"
# IMAGE_PATH        = "../data/spatial/tissue_hires_image.png"
# SCALEFACTOR_PATH  = "../data/binned_outputs/square_008um/spatial/scalefactors_json.json"
# POINT_SIZE        = 0.3
#
# # ── 1. Capture-area overlay on H&E ──
# df  = pd.read_csv(CSV_PATH)
# img = plt.imread(IMAGE_PATH)
# with open(SCALEFACTOR_PATH) as f:
#     scale_factors = json.load(f)
#
# hires_scale = scale_factors["tissue_hires_scalef"]
# scaled_x = df["X_Coordinate"] * hires_scale
# scaled_y = df["Y_Coordinate"] * hires_scale
#
# min_x, max_x = scaled_x.min(), scaled_x.max()
# min_y, max_y = scaled_y.min(), scaled_y.max()
#
# fig, ax = plt.subplots(figsize=(12, 12))
# ax.imshow(img)
# rect = patches.Rectangle(
#     (min_x, min_y), max_x - min_x, max_y - min_y,
#     linewidth=3, edgecolor="lime", facecolor="none",
#     linestyle="--", label="Sequenced Capture Area (6.5×6.5 mm)")
# ax.add_patch(rect)
# ax.scatter(scaled_x, scaled_y, color="cyan", alpha=0.05, s=1, edgecolors="none")
# ax.set_title("Mapping the Visium HD Capture Area to the Full Tissue Scan", fontsize=16)
# ax.legend(loc="upper right", fontsize=12)
# ax.axis("off")
# plt.tight_layout()
# plt.show()
#
# # ── 2. Side-by-side: full tissue vs. sketch ──
# all_x    = df["X_Coordinate"].values
# all_y    = df["Y_Coordinate"].values
# sketched = df["Is_Sketched"].values.astype(bool)
# n_total  = len(df)
# n_sketch = int(sketched.sum())
# print(f"Loaded {n_total} bins, {n_sketch} sketched ({100*n_sketch/n_total:.2f}%)")
#
# fig2, (ax_full, ax_sketch) = plt.subplots(1, 2, figsize=(16, 7))
#
# ax_full.scatter(all_x, all_y, s=POINT_SIZE, c="steelblue", alpha=0.6)
# ax_full.set_title(f"Full Tissue  (N = {n_total:,})", fontsize=13)
# ax_full.set_xlabel("X"); ax_full.set_ylabel("Y")
# ax_full.set_aspect("equal"); ax_full.invert_yaxis()
#
# ax_sketch.scatter(all_x[sketched], all_y[sketched], s=POINT_SIZE, c="crimson", alpha=0.8)
# ax_sketch.set_title(f"Sketch  (k = {n_sketch:,},  {100*n_sketch/n_total:.2f}%)", fontsize=13)
# ax_sketch.set_xlabel("X"); ax_sketch.set_ylabel("Y")
# ax_sketch.set_aspect("equal"); ax_sketch.invert_yaxis()
#
# for ax in (ax_full, ax_sketch):
#     ax.set_xlim(all_x.min() - 10, all_x.max() + 10)
#     ax.set_ylim(all_y.max() + 10, all_y.min() - 10)
#
# fig2.suptitle("Spatially-Aware Sketching: Full Tissue vs. Downsampled Sketch",
#               fontsize=15, fontweight="bold")
# fig2.tight_layout(rect=[0, 0, 1, 0.95])
# plt.show()
#
# ══════════════════════════════════════════════════════════════════════

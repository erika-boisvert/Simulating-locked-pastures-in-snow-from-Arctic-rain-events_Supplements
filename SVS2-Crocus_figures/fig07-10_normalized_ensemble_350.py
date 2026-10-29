import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import pandas as pd
import matplotlib as mpl
import sys
import time

mpl.rcParams.update({
    # Font family: 'serif', 'sans-serif', 'monospace', 'cursive', 'fantasy'
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Times New Roman', 'Helvetica', 'DejaVu Sans'],  # fallback order
    'font.size':         10,      # default for all text
    'font.style':        'normal', # 'normal', 'italic', 'oblique'
    'font.weight':       'normal', # 'normal', 'bold', 'light', 400, 700, etc.

    # Override specific elements independently
    'axes.titlesize':    14,
    'axes.titleweight':  'bold',
    'axes.labelsize':    10,
    'axes.labelweight':  'normal',

    'xtick.labelsize':   10,
    'ytick.labelsize':   10,

    'legend.fontsize':   10,
    'legend.title_fontsize': 11,

    'figure.titlesize':  14,      # suptitle
    'figure.titleweight':'bold',
    })

if len(sys.argv) != 3:
    raise ValueError('Requires two argument: <ensemble folder> <year>')

# Land surface scheme : select among svs1 and svs2
ens_name = sys.argv[-2]
year = sys.argv[-1]

# Define parameters
result_folder = os.path.join(r'/home/ebv000/experiments', ens_name)
output_folder = os.path.join(result_folder, 'all_outputs')
# os.makedirs(output_folder, exist_ok=True)
output_filename = os.path.join(output_folder, 'ensemble_snowpack_dens350_normalized' + str(year) + '.nc')

start_date = str(year) + '-07-02'
end_date = str(int(year)+1) + '-06-30'
density_threshold = 350

# Resolution of the normalized (0-1) depth axis. 1000 bins gives 0.1%-of-depth
# resolution, comparable to the ~1mm resolution the previous absolute-depth
# version used for a ~0.3-0.5 m snowpack.
N_DEPTH_BINS = 1000


# ------------------------------------------------------------------------- #
# --- Find all NetCDF files for this ensemble ---
#
# The per-member archives now live in a single consolidated folder named after
# the simulation (rather than scattered loose folders). We don't assume a
# single fixed layout: we search a short list of the most likely locations,
# recursively, and report exactly what we found so a wrong guess is obvious
# immediately rather than silently returning zero files.
# ------------------------------------------------------------------------- #
CANDIDATE_BASE_DIRS = [
    os.path.join(r'/home/ebv000/work', ens_name + '_outputs'),           # consolidated ensemble folder (most likely)
    # os.path.join(r'/home/ebv000/experiments', ens_name, 'all_outputs'),  # original assumed layout
]

nc_files = []
found_base = None
for base_dir in CANDIDATE_BASE_DIRS:
    if not os.path.isdir(base_dir):
        print(f"  (skipping, does not exist: {base_dir})")
        continue
    candidates = sorted(glob.glob(os.path.join(base_dir, '*.nc'), recursive=True))
    # Exclude any previously-generated ensemble output file sitting in the same tree
    candidates = [f for f in candidates if 'ensemble_snowpack_dens' not in os.path.basename(f)]
    print(f"  Searched {base_dir} -> {len(candidates)} file(s)")
    if candidates:
        nc_files = candidates
        found_base = base_dir
        break

if not nc_files:
    print("No NetCDF files found in any candidate location:")
    for d in CANDIDATE_BASE_DIRS:
        print(f"  - {d}")
    print("Update CANDIDATE_BASE_DIRS at the top of this script if the ensemble "
          "was archived somewhere else.")
    sys.exit(1)

print(f"Found {len(nc_files)} NetCDF files under {found_base}:")
for f in nc_files:
    print(f"  - {f}")


# ------------------------------------------------------------------------- #
def centers_to_edges(points):
    """Turn an array of bin *centers* into an array of bin *edges* (len+1),
    using the midpoint between neighbours and linear extrapolation at the ends.
    Shared by the time axis and (in the plotting function) the mean-depth axis.
    """
    points = np.asarray(points, dtype=np.float64)
    edges = np.zeros(len(points) + 1)
    if len(points) > 1:
        edges[0] = points[0] - (points[1] - points[0]) / 2
        edges[1:-1] = (points[:-1] + points[1:]) / 2
        edges[-1] = points[-1] + (points[-1] - points[-2]) / 2
    elif len(points) == 1:
        edges[0] = points[0] - 0.5
        edges[-1] = points[0] + 0.5
    return edges


def simulation_name_from_path(file_path, base_dir):
    """Derive a readable run label from a file path, relative to the base
    directory it was found under (e.g. 'member03/out_svs2.nc' rather than
    the full absolute path). Falls back to the bare filename if base_dir
    doesn't prefix the path for some reason.
    """
    try:
        rel = os.path.relpath(file_path, base_dir)
    except ValueError:
        rel = os.path.basename(file_path)
    return rel


def create_normalized_density_netcdf(file_paths, output_filename, start_date='2003-07-02',
                                      end_date='2004-06-30', density_threshold=350,
                                      n_depth_bins=N_DEPTH_BINS, base_dir=None):
    """
    Create a NetCDF of snow density packed onto a *normalized* depth axis.

    Instead of an absolute-depth (mm) grid, each run's layers are expressed as
    a fraction of THAT run's own total snow depth at that timestep (0 = ground,
    1 = surface). This puts every run on the same 0-1 scale regardless of its
    actual depth, so the internal density banding can be compared as ratios
    across runs. `total_snow_depth` (absolute, in meters) is still stored per
    run/time so the plotting step can redistribute the normalized profile onto
    real depth units using the ensemble-average depth.

    Optimizations relative to the original absolute-depth version are kept:
      - Each input file is opened exactly once.
      - Per-timestep/per-layer values are pulled out as plain numpy arrays up
        front instead of repeated label-based indexing inside nested loops.
      - `thickness <= 0` marks unused layers (this dataset fills them with
        mass=0, not NaN), so processing stops at the real snow surface.
      - Layer boundaries within a timestep use a vectorized `cumsum`.
    """
    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)

    print("Step 1: Loading all files (single pass)...")
    t0 = time.time()

    all_times = None
    per_run_data = []  # list of dicts: {thickness, density, depth} numpy arrays, per run

    for file_path in file_paths:
        try:
            ds = xr.open_dataset(file_path)
            ds_time = ds.sel(time=slice(start_dt, end_dt))
            if len(ds_time.time) == 0:
                ds.close()
                per_run_data.append(None)
                continue

            if all_times is None:
                all_times = ds_time.time.values

            thickness = (ds_time['SNOMA_ML'] / ds_time['SNODEN_ML']).values.astype(np.float32)
            density = ds_time['SNODEN_ML'].values.astype(np.float32)
            depth = (ds_time['SNODP'].values.astype(np.float32)
                     if 'SNODP' in ds.variables else np.zeros(len(ds_time.time), dtype=np.float32))

            # Treat non-positive thickness as "no snow in this layer" (this dataset
            # uses 0, not NaN, for unused layers).
            thickness = np.where(thickness > 0, thickness, 0.0)

            per_run_data.append({'thickness': thickness, 'density': density, 'depth': depth})
            ds.close()
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            per_run_data.append(None)
            continue

    if all_times is None:
        print("No valid data found in the specified time range")
        return None

    n_times = len(all_times)
    n_runs = len(file_paths)
    print(f"  Loaded {n_runs} run(s) in {time.time() - t0:.1f}s")
    print(f"Time steps: {n_times} hourly steps")

    # Normalized depth axis: index i represents the fraction i/n_depth_bins of
    # each run's own total depth at that timestep (0 = ground, 1 = surface).
    depth_fraction = np.arange(0, n_depth_bins + 1) / n_depth_bins
    n_depth_pts = len(depth_fraction)
    print(f"Depth resolution: {n_depth_bins} bins (fraction 0-1 of each run's own total depth)")

    output_density = np.full((n_runs, n_times, n_depth_pts), np.nan, dtype=np.float32)
    total_snow_depth_array = np.full((n_runs, n_times), np.nan, dtype=np.float32)

    print("Step 2: Packing each run onto the normalized (0-1) depth grid...")
    t0 = time.time()

    for run_idx, data in enumerate(per_run_data):
        if data is None:
            continue
        thickness_all = data['thickness']   # (n_times, n_layers)
        density_all = data['density']       # (n_times, n_layers)
        depth_all = data['depth']           # (n_times,)
        total_snow_depth_array[run_idx, :] = depth_all

        for t_idx in range(n_times):
            thickness_row = thickness_all[t_idx]
            valid = thickness_row > 0
            if not valid.any():
                continue

            total_depth = depth_all[t_idx]
            if np.isnan(total_depth) or total_depth <= 0:
                # No meaningful total depth to normalize against -- skip.
                continue

            valid_thickness = thickness_row[valid]
            valid_density = density_all[t_idx][valid]

            # Vectorized cumulative boundaries: layer 0 = surface (top), last
            # valid layer = ground. cum[0] = 0, cum[k] = sum of first k thicknesses.
            cum = np.concatenate(([0.0], np.cumsum(valid_thickness)))
            height_top = total_depth - cum[:-1]
            height_bottom = total_depth - cum[1:]

            # Normalize to a fraction of THIS run's own total depth at this
            # time -- this is what puts every run on the same 0-1 scale
            # regardless of its actual absolute depth.
            frac_top = np.clip(height_top / total_depth, 0.0, 1.0)
            frac_bottom = np.clip(height_bottom / total_depth, 0.0, 1.0)

            bottom_idx = np.clip(np.round(frac_bottom * n_depth_bins).astype(int), 0, n_depth_pts - 1)
            top_idx = np.clip(np.round(frac_top * n_depth_bins).astype(int), 0, n_depth_pts - 1)

            row = np.full(n_depth_pts, np.nan, dtype=np.float32)
            for b, t, d in zip(bottom_idx, top_idx, valid_density):
                if b < t:
                    row[b:t] = d
            output_density[run_idx, t_idx, :] = row

        if (run_idx + 1) % max(1, n_runs // 10) == 0 or run_idx == n_runs - 1:
            print(f"  Processed run {run_idx + 1}/{n_runs}")

    print(f"  Packing complete in {time.time() - t0:.1f}s")

    print("Step 3: Creating output NetCDF...")
    run_names = np.array([simulation_name_from_path(fp, base_dir or '') for fp in file_paths])

    output_ds = xr.Dataset({
        'snow_density': (['model_run', 'time', 'depth_fraction'], output_density)
    }, coords={
        'model_run': np.arange(1, n_runs + 1),
        'time': all_times,
        'depth_fraction': depth_fraction,
        'simulation': ('model_run', run_names),
    })
    output_ds['snow_density'].attrs = {
        'long_name': 'Snow density on a normalized 0-1 depth axis '
                     '(0 = ground, 1 = surface, relative to each run\'s own total depth)',
        'units': 'kg/m3'
    }
    output_ds.to_netcdf(output_filename)
    print(f"Saved density data to {output_filename}")

    print("Step 4: Creating threshold counts...")
    high_density_mask = output_density > density_threshold
    count_matrix = np.sum(high_density_mask, axis=0).astype(np.float32)

    # Cells where NO run has any snow at all (every run is NaN here) are genuine
    # "no data" -- keep these as NaN so the plot leaves them transparent (showing
    # the ggplot background) rather than opaque white, which is reserved for cells
    # where at least one run *does* have snow but none exceed the threshold.
    any_data = np.any(~np.isnan(output_density), axis=0)  # (n_times, n_depth_pts)
    count_matrix[~any_data] = np.nan

    count_ds = xr.Dataset({
        'high_density_count': (['time', 'depth_fraction'], count_matrix),
        'total_snow_depth': (['model_run', 'time'], total_snow_depth_array)
    }, coords={
        'time': all_times,
        'depth_fraction': depth_fraction,
        'model_run': np.arange(1, n_runs + 1),
        'simulation': ('model_run', run_names),
    })
    count_ds['high_density_count'].attrs = {
        'long_name': f'Number of model runs with snow density > {density_threshold} kg/m3 '
                     'at this normalized depth fraction (NaN = no run has snow here)',
        'units': 'count'
    }
    count_ds['total_snow_depth'].attrs = {
        'long_name': 'Total (absolute) snow depth for each model run -- kept so the '
                     'normalized profile can later be redistributed onto real depth units',
        'units': 'm'
    }
    count_filename = output_filename.replace('.nc', '_counts.nc')
    count_ds.to_netcdf(count_filename)
    print(f"Saved count data to {count_filename}")

    return output_ds, count_ds


# ------------------------------------------------------------------------- #
def plot_from_netcdf(count_filename, start_date='2003-07-02', end_date='2004-06-30',
                      density_threshold=350):
    """
    Plot the normalized density-count grid, redistributed onto real depth units
    using the ensemble-average total snow depth at each time. Every run's
    profile is on the same 0-1 fraction scale internally (so the banding is
    directly comparable across runs), but the vertical axis you see is scaled
    by mean_depth(t) = average of `total_snow_depth` across all runs at that
    timestep, so the plot still reads in physically meaningful meters and the
    overall envelope traces the ensemble-average snow depth over the season.
    """
    print("Loading pre-processed data...")
    count_ds = xr.open_dataset(count_filename)

    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)
    count_data = count_ds.sel(time=slice(start_dt, end_dt))

    times = count_data.time.values
    depth_fraction = count_data.depth_fraction.values

    count_matrix = np.ma.masked_invalid(count_data['high_density_count'].values.T)  # (depth_fraction, time)
    total_snow_depth_all = count_data['total_snow_depth'].values  # (model_run, time), absolute meters

    # The ensemble-average total depth at each time -- this is what the
    # normalized profile gets redistributed onto for display.
    mean_depth_per_time = np.nanmean(total_snow_depth_all, axis=0)  # (n_times,)

    plt.style.use('ggplot')

    fig, ax = plt.subplots(1, 1, figsize=(11, 6))

    time_centers = mpl.dates.date2num(times)
    time_edges = centers_to_edges(time_centers)

    # Normalized-depth bin edges (left-edge style, matching the 0..1 sample points).
    frac_edges = np.zeros(len(depth_fraction) + 1)
    frac_edges[:-1] = depth_fraction
    if len(depth_fraction) > 1:
        frac_edges[-1] = depth_fraction[-1] + (depth_fraction[1] - depth_fraction[0])
    else:
        frac_edges[-1] = depth_fraction[-1] + 0.001

    # Ensemble-mean-depth edges (center-style, matching how `time_edges` was built
    # from `time_centers`) -- never let this go negative.
    mean_depth_edges = np.clip(centers_to_edges(mean_depth_per_time), 0, None)

    # Build the "fan" mesh: each time-column is scaled by the ensemble-average
    # depth at that time, so DEPTH2D[i, j] = frac_edges[i] * mean_depth_edges[j].
    TIME2D = np.tile(time_edges[None, :], (len(frac_edges), 1))
    DEPTH2D = frac_edges[:, None] * mean_depth_edges[None, :]

    n_runs = total_snow_depth_all.shape[0]

    cmap = plt.cm.CMRmap_r
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = (1.0, 1.0, 1.0, 1.0)  # zero-count cells: opaque white (was transparent)
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

    # bounds = np.linspace(0, n_runs, n_runs + 1)
    bounds = np.linspace(0, 100, n_runs + 1)
    bounds_ticks = np.linspace(0, 100, 11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    if str(year) == '2003':
        highlight_start1 = pd.to_datetime('2003-10-02T12:00:00')
        highlight_end1 = pd.to_datetime('2003-10-03T00:00:00')
        highlight_start2 = pd.to_datetime('2003-10-04T00:00:00')
        highlight_end2 = pd.to_datetime('2003-10-05T12:00:00')

        highlight_x1 = mpl.dates.date2num(highlight_start1)
        highlight_x2 = mpl.dates.date2num(highlight_end1)
        highlight_x3 = mpl.dates.date2num(highlight_start2)
        highlight_x4 = mpl.dates.date2num(highlight_end2)
        ax.axvspan(highlight_x1, highlight_x2, color='red', alpha=0.16, zorder=1)
        ax.axvspan(highlight_x3, highlight_x4, color='red', alpha=0.16, zorder=1)

    ax.pcolormesh(TIME2D, DEPTH2D, count_matrix, cmap=cmap, norm=norm, shading='flat', zorder=2)


    # Background silhouette = the ensemble-average depth curve itself, since
    # that's exactly what the color mesh's top boundary follows.
    ax.fill_between(time_centers, mean_depth_per_time, 0, color='whitesmoke', zorder=1.5)
    # for i in range(n_runs):
    #     ax.plot(time_centers, total_snow_depth_all[i, :],
    #             color='steelblue', alpha=0.9, linewidth=1, linestyle='-.', zorder=3)
    ax.plot(time_centers, mean_depth_per_time, color='steelblue', linewidth=1.5, linestyle='--',
            zorder=3.5, label='Ensemble-mean depth')

    # --- Inset ---
    # Anchored to the main axes (not the figure) so it aligns with the y-axis (left
    # spine) and the axes' own top border -- which sits below the title -- rather
    # than the figure's outer edge.
    INSET_WIDTH, INSET_HEIGHT = 0.30, 0.65  # axes-fraction; height controls how low the bottom edge sits
    ax_inset = ax.inset_axes([0.0, 1.0 - INSET_HEIGHT, INSET_WIDTH, INSET_HEIGHT])

    zoom_start = pd.to_datetime(str(year) + '-09-29')
    zoom_end = pd.to_datetime(str(year) + '-10-13')

    zoom_mask = (times >= np.datetime64(zoom_start)) & (times <= np.datetime64(zoom_end))
    zoom_time_centers = time_centers[zoom_mask]
    zoom_mean_depth = mean_depth_per_time[zoom_mask]
    zoom_count_matrix = count_matrix[:, zoom_mask]

    zoom_time_edges = centers_to_edges(zoom_time_centers)
    zoom_mean_depth_edges = np.clip(centers_to_edges(zoom_mean_depth), 0, None)

    ZOOM_TIME2D = np.tile(zoom_time_edges[None, :], (len(frac_edges), 1))
    ZOOM_DEPTH2D = frac_edges[:, None] * zoom_mean_depth_edges[None, :]

    ax_inset.pcolormesh(ZOOM_TIME2D, ZOOM_DEPTH2D, zoom_count_matrix,
                         cmap=cmap, norm=norm, shading='flat', zorder=2)

    zoom_x1_num = mpl.dates.date2num(zoom_start)
    zoom_x2_num = mpl.dates.date2num(zoom_end)
    
    if str(year) == '2003':
        if highlight_x1 < zoom_x2_num and highlight_x2 > zoom_x1_num:
            ax_inset.axvspan(max(highlight_x1, zoom_x1_num), min(highlight_x2, zoom_x2_num),
                            color='red', alpha=0.16, zorder=1)
        if highlight_x3 < zoom_x2_num and highlight_x4 > zoom_x1_num:
            ax_inset.axvspan(max(highlight_x3, zoom_x1_num), min(highlight_x4, zoom_x2_num),
                            color='red', alpha=0.16, zorder=1)

    ax_inset.fill_between(zoom_time_centers, zoom_mean_depth, 0, color='whitesmoke', zorder=1.5)
    # for i in range(n_runs):
    #     ax_inset.plot(zoom_time_centers, total_snow_depth_all[i, zoom_mask],
    #                   color='steelblue', alpha=0.9, linewidth=1, linestyle='-.', zorder=3)
    ax_inset.plot(zoom_time_centers, zoom_mean_depth, color='steelblue', linewidth=1.5,
                  linestyle='--', zorder=3.5)

    ax_inset.set_xlim(zoom_time_edges[0], zoom_time_edges[-1])
    # A little headroom above whichever is taller in this window: the mean-depth
    # envelope or an individual run's curve.
    zoom_ylim_max = np.nanmax([np.nanmax(zoom_mean_depth_edges),
                                np.nanmax(total_snow_depth_all[:, zoom_mask])])
    zoom_ylim_max = zoom_ylim_max * 1.1 if np.isfinite(zoom_ylim_max) and zoom_ylim_max > 0 else 0.1
    ax_inset.set_ylim(0, zoom_ylim_max)
    ax_inset.tick_params(axis='both', which='both', bottom=False, left=False,
                          labelbottom=False, labelleft=False)
    ax_inset.grid(True, zorder=1)
    for spine in ax_inset.spines.values():
        spine.set_color('black')

    # Rectangle on the main plot + connector lines to the inset are both derived
    # automatically from ax_inset's actual position and its xlim/ylim -- no manual
    # coordinates to keep in sync if the inset ever moves.
    ax.indicate_inset_zoom(ax_inset, edgecolor='k', linewidth=1, alpha=0.8, zorder=4)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Average Snow Depth (m)', fontsize=12)
    # ax.set_title(f'{ens_name}: Normalized Model Runs Redistributed onto Ensemble-Average Depth\n'
    #              f'(Density > {density_threshold} kg/m³, {start_date} to {end_date})',
    #              fontsize=14)
    # ax.legend(loc='upper right', frameon=True, fontsize=9)

    ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b %Y'))  # e.g. "Jul 2003"
    ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(interval=1))
    ax.tick_params(axis='x', which='major', length=7, width=1.5, color='black')
    ax.tick_params(axis='x', which='minor', length=3, width=0.8, color='black')

    # Y-limit derived from the data itself (ensemble-mean envelope and the
    # individual run curves) rather than a hardcoded value tuned to one dataset.
    main_ylim_max = np.nanmax([np.nanmax(total_snow_depth_all), np.nanmax(mean_depth_edges)])
    main_ylim_max = main_ylim_max * 1.05 if np.isfinite(main_ylim_max) and main_ylim_max > 0 else 0.5
    ax.set_ylim(0, main_ylim_max)

    # Frame around just the data area (ggplot style otherwise hides/lightens some
    # spines) -- doesn't touch the title, colorbar, or tick labels, since those all
    # live outside this box.
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1)

    ax2 = fig.add_axes([1.015, 0.1, 0.02, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
                                    spacing='proportional', ticks=bounds_ticks, boundaries=bounds, format='%1i')
    ax2.set_ylabel('Probability of high density snow (> 350 kg m$^{-3}$)', size=12)

    # Frame around the colorbar/legend, matching the main plot's frame.
    for spine in ax2.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1)

    ax.grid(True, zorder=1)
    plt.xticks(rotation=45)
    plt.tight_layout()

    count_ds.close()
    return fig, ax


# ------------------------------------------------------------------------- #
if os.path.exists(output_filename):
    print(f"Output file {output_filename} already exists.")
    use_existing = input("Use existing file? (y/n): ").lower().strip()
    if use_existing == 'n':
        create_normalized_density_netcdf(nc_files, output_filename, start_date, end_date,
                                          density_threshold, base_dir=found_base)
else:
    create_normalized_density_netcdf(nc_files, output_filename, start_date, end_date,
                                      density_threshold, base_dir=found_base)

count_filename = output_filename.replace('.nc', '_counts.nc')
fig, ax = plot_from_netcdf(count_filename, start_date, end_date, density_threshold)
plt.savefig(os.path.join(result_folder, 'snowpack_density_normalized_' + str(year) + '.png'), dpi=300, bbox_inches='tight')
plt.show()

print("Processing complete!")
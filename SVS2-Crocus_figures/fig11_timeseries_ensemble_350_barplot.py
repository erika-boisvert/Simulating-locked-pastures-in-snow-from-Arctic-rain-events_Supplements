"""
Ensemble version of Figure 10: snowpack density on 31 October of every year,
shown as the probability (% of ensemble members) that a layer is denser than
350 kg m-3.

One figure is produced (REDISTRIBUTED): each member's profile is normalized
(0 = ground, 1 = snow surface), then stretched back onto real depth using the
ensemble-mean snow depth of that year on 31 October. The y-axis is therefore
in metres. (The normalized-axis version is still available by calling
plot_ensemble(ds, redistribute=False, outfile=...).)

Additional outputs (CSV):
  - basal_20mm_probabilities.csv : probability of density > threshold in the
                                   bottom 20 mm, per year
  - snow_onset_ensemble_mean.csv : ensemble-mean snow onset date, per year
  - snow_onset_all_members.csv   : onset date of every member and year

Snow onset = start of the FINAL persistent snow cover of the winter, i.e. the
date from which snow stays on the ground (apart from short gaps) until spring
melt-out. If the snowpack melts out completely in November and then
re-accumulates in December and persists until June, the December date is
reported.

Usage:
    python normalized_ensemble_350_oct31.py <ensemble name> [--recompute] [--no-onset]

<ensemble name> is the sim_name used in Crocus_ensemble.sh for the long
(1980-2020) ensemble. Member files are read from
    /home/ebv000/work/<ensemble name>_outputs/*.nc
exactly like normalized_ensemble_350.py (Figure 7).
"""
import argparse
import glob
import os
import sys
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D
from matplotlib.path import Path
from matplotlib.patches import Rectangle

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Times New Roman', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
})

# ------------------------------------------------------------------------- #
# Settings
# ------------------------------------------------------------------------- #
FIRST_YEAR = 1985          # first 5 years (1980-84) are spin-up
LAST_YEAR = 2020           # NB: 31 Oct 2020 must exist in the simulation
SNAPSHOT_MONTH_DAY = '10-31'
SNAPSHOT_HOUR = 12         # which hour of 31 Oct to sample (nearest time step is used)
MAX_TIME_OFFSET = pd.Timedelta('1D')   # error out (skip year) if no step this close

DENSITY_THRESHOLD = 350    # kg m-3
N_DEPTH_BINS = 1000        # resolution of the normalized (0-1) depth axis
COLUMN_WIDTH = 1.0         # width of each year's column (1 = columns touch)
SIDE_MARGIN = 0.75         # empty space left/right of the bars, in bar widths
TICK_STEP_YEARS = 3        # major (labelled) tick every N bars; minor tick on every bar
OUTLINE_COLOR = 'steelblue'

# Basal layer analysis
BASAL_THICKNESS_MM = 20
BASAL_YEARS_OF_INTEREST = [1988, 1996, 1998, 2003, 2008]

# Snow onset definition: start of the FINAL persistent snow-covered period of
# the winter season. The season runs from ONSET_SEARCH_START of year Y to
# ONSET_SEARCH_END of year Y+1. Snow-free interruptions shorter than
# ONSET_GAP_DAYS are ignored (brief melt / wind-scour dropouts), and a
# snow-covered period must last at least ONSET_MIN_SEASON_DAYS to count as
# "the winter snowpack" (so a late-summer dusting is not picked up).
ONSET_DEPTH_THRESHOLD = 0.01     # m; snow depth must exceed this to count as snow cover
ONSET_GAP_DAYS = 3               # snow-free spells >= this length end a snow period
ONSET_MIN_SEASON_DAYS = 30       # minimum length of the persistent winter snow period
ONSET_SEARCH_START = '08-01'     # season start (month-day, year Y)
ONSET_SEARCH_END = '07-31'       # season end (month-day, year Y+1)
ONSET_PROVISIONAL_BEFORE = '05-31'   # if a run ends before this date in Y+1, flag the year

# Markers read from the published Figure 10 (droplets = positive wet-snow
# detections at Thomsen River; star = ROS-associated ungulate die-off, drawn
# on top of the droplet stack). Replace with the exact counts from
# Langlois et al. (2026) if they differ.
#   1987: star only          1988: 1 droplet + star    1990: star only
#   1993: star only          1996: 1 droplet           1998: 1 droplet
#   2003: 4 droplets + star  2008: 1 droplet
WET_SNOW_DETECTIONS = {1988: 1, 1996: 1, 1998: 1, 2003: 4, 2008: 1}   # year: n detections
ROS_DIEOFF_YEARS = [1987, 1988, 1990, 1993, 2003]

VARS = ['SNOMA_ML', 'SNODEN_ML', 'SNODP']

# ------------------------------------------------------------------------- #
# Arguments / paths
# ------------------------------------------------------------------------- #
parser = argparse.ArgumentParser()
parser.add_argument('ens_name', help='ensemble (sim_name) folder')
parser.add_argument('--recompute', action='store_true',
                    help='rebuild the intermediate NetCDF even if it exists')
parser.add_argument('--no-onset', action='store_true',
                    help='skip the snow onset table (it re-reads all member files)')
args = parser.parse_args()
ens_name = args.ens_name

result_folder = os.path.join('/home/ebv000/experiments', ens_name)
output_folder = os.path.join(result_folder, 'all_outputs')
data_filename = os.path.join(output_folder, 'ensemble_snowpack_dens350_oct31_normalized.nc')
CANDIDATE_BASE_DIRS = [os.path.join('/home/ebv000/work', ens_name + '_outputs')]


def find_member_files():
    for base in CANDIDATE_BASE_DIRS:
        if not os.path.isdir(base):
            print(f"  (skipping, does not exist: {base})")
            continue
        files = sorted(glob.glob(os.path.join(base, '*.nc')))
        files = [f for f in files if 'ensemble_snowpack_dens' not in os.path.basename(f)]
        print(f"  Searched {base} -> {len(files)} file(s)")
        if files:
            return files, base
    sys.exit("No member NetCDF files found. Edit CANDIDATE_BASE_DIRS.")


# ------------------------------------------------------------------------- #
# Step 1: extract 31 Oct profiles, pack on normalized depth axis
# ------------------------------------------------------------------------- #
def build_dataset(nc_files, base_dir):
    years = np.arange(FIRST_YEAR, LAST_YEAR + 1)
    targets = pd.to_datetime([f'{y}-{SNAPSHOT_MONTH_DAY} {SNAPSHOT_HOUR:02d}:00' for y in years])
    n_runs, n_years, N = len(nc_files), len(years), N_DEPTH_BINS

    density_norm = np.full((n_runs, n_years, N), np.nan, dtype=np.float32)
    total_depth = np.full((n_runs, n_years), np.nan, dtype=np.float32)
    run_names = []

    t0 = time.time()
    for r, fp in enumerate(nc_files):
        run_names.append(os.path.basename(fp).replace('.nc', ''))
        try:
            with xr.open_dataset(fp) as ds:
                # Pull ONLY the 31-Oct time steps out of the 40-year hourly file
                sub = ds[VARS].sel(time=targets, method='nearest').load()
        except Exception as e:
            print(f"  Error reading {fp}: {e}")
            continue

        offset = np.abs(sub.time.values - targets.values)
        n_t = len(years)
        thick_all = (sub['SNOMA_ML'] / sub['SNODEN_ML']).values.reshape(n_t, -1)
        dens_all = sub['SNODEN_ML'].values.reshape(n_t, -1)
        depth_all = sub['SNODP'].values.reshape(n_t)

        for j in range(n_t):
            if offset[j] > MAX_TIME_OFFSET.to_timedelta64():
                continue                                  # year not covered by this run
            D = depth_all[j]
            th = np.nan_to_num(thick_all[j], nan=0.0)
            th = np.where(th > 0, th, 0.0)
            valid = th > 0
            total_depth[r, j] = D
            if not valid.any() or not np.isfinite(D) or D <= 0:
                continue                                  # no snow on that date

            cum = np.concatenate(([0.0], np.cumsum(th[valid])))
            top = np.clip((D - cum[:-1]) / D, 0, 1)       # layer top, fraction of depth
            bot = np.clip((D - cum[1:]) / D, 0, 1)        # layer bottom
            b_idx = np.round(bot * N).astype(int)
            t_idx = np.round(top * N).astype(int)
            row = np.full(N, np.nan, dtype=np.float32)
            for b, t, d in zip(b_idx, t_idx, dens_all[j][valid]):
                if b < t:
                    row[b:t] = d
            density_norm[r, j] = row

        if (r + 1) % max(1, n_runs // 10) == 0 or r == n_runs - 1:
            print(f"  Processed {r + 1}/{n_runs} runs ({time.time() - t0:.0f}s)")

    # Probability of density > threshold at each normalized height, per year.
    # Denominator = members that actually have snow in that year on that date.
    valid_cnt = np.sum(np.isfinite(density_norm), axis=0)                  # (year, bin)
    high_cnt = np.sum(density_norm > DENSITY_THRESHOLD, axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        prob = np.where(valid_cnt > 0, 100.0 * high_cnt / valid_cnt, np.nan)
    n_with_snow = np.sum(np.isfinite(density_norm[:, :, 0]), axis=0)       # (year,)

    ds_out = xr.Dataset(
        {
            'prob_high_density': (['year', 'depth_fraction'], prob.astype(np.float32),
                                  {'units': '%', 'long_name': f'% of members with density > {DENSITY_THRESHOLD} kg/m3'}),
            'n_members_with_snow': (['year'], n_with_snow),
            'total_snow_depth': (['model_run', 'year'], total_depth, {'units': 'm'}),
            'snow_density': (['model_run', 'year', 'depth_fraction'], density_norm, {'units': 'kg/m3'}),
        },
        coords={'year': years,
                'depth_fraction': (np.arange(N) + 0.5) / N,     # bin centres
                'model_run': np.arange(1, n_runs + 1),
                'simulation': ('model_run', np.array(run_names))},
    )
    ds_out.to_netcdf(data_filename)
    print(f"Saved {data_filename}")
    return ds_out


# ------------------------------------------------------------------------- #
# Basal-layer values (as drawn in the redistributed figure)
# ------------------------------------------------------------------------- #
def print_basal_values(ds, thickness_mm=20, years_of_interest=None, csv_path=None):
    """Probability (%) of density > DENSITY_THRESHOLD in the bottom
    `thickness_mm` of the snowpack, using the same ensemble-mean depth
    stretching as the redistributed figure."""
    years = ds.year.values
    prob = ds['prob_high_density'].values                            # (year, bin)
    frac_centres = ds['depth_fraction'].values
    mean_depth = np.nanmean(ds['total_snow_depth'].values, axis=0)   # (year,) in m
    h = thickness_mm / 1000.0

    print(f"\nBottom {thickness_mm} mm: P(density > {DENSITY_THRESHOLD} kg/m3) [%], "
          f"redistributed onto ensemble-mean depth")
    print(f"{'Year':>5} {'mean depth (mm)':>16} {'layer used (mm)':>16} "
          f"{'mean %':>8} {'min %':>7} {'max %':>7} {'base bin %':>11}")

    rows = []
    for j, yr in enumerate(years):
        D = mean_depth[j]
        if not np.isfinite(D) or D <= 0 or not np.isfinite(prob[j]).any():
            print(f"{yr:>5} {'n/a':>16}")
            rows.append((yr, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))
            continue
        frac = min(h / D, 1.0)                 # if the pack is thinner than 20 mm, use all of it
        sel = frac_centres < frac
        if not sel.any():
            sel[0] = True
        vals = prob[j, sel]
        mean_v, min_v, max_v = np.nanmean(vals), np.nanmin(vals), np.nanmax(vals)
        base_v = prob[j, 0]
        used_mm = frac * D * 1000
        rows.append((yr, D * 1000, used_mm, mean_v, min_v, max_v, base_v))
        flag = '  <--' if years_of_interest and yr in years_of_interest else ''
        print(f"{yr:>5} {D * 1000:>16.1f} {used_mm:>16.1f} "
              f"{mean_v:>8.1f} {min_v:>7.1f} {max_v:>7.1f} {base_v:>11.1f}{flag}")

    if csv_path:
        pd.DataFrame(rows, columns=['year', 'mean_depth_mm', 'layer_used_mm',
                                    'prob_mean_pct', 'prob_min_pct', 'prob_max_pct',
                                    'prob_base_bin_pct']).to_csv(csv_path, index=False)
        print(f"Saved {csv_path}")


# ------------------------------------------------------------------------- #
# Snow onset (mean over all ensemble members)
# ------------------------------------------------------------------------- #
def _final_snow_period_start(snow, index, gap_n, min_n):
    """Start time of the LAST snow-covered period lasting >= min_n steps.
    `snow` is a boolean array (True = snow on the ground). Snow-free spells
    shorter than gap_n steps do not split a period; spells of gap_n steps or
    more do. Returns None if no qualifying period exists."""
    idx = np.flatnonzero(snow)
    if idx.size == 0:
        return None
    breaks = np.flatnonzero(np.diff(idx) - 1 >= gap_n)     # long snow-free gaps
    starts = np.r_[idx[0], idx[breaks + 1]]                # first snow step of each period
    ends = np.r_[idx[breaks], idx[-1]]                     # last snow step of each period
    ok = np.flatnonzero((ends - starts + 1) >= min_n)      # periods long enough to count
    if ok.size == 0:
        return None
    return index[starts[ok[-1]]]                           # the last qualifying period


def compute_snow_onset_table(nc_files, csv_path, members_csv_path=None):
    """Snow onset per member and year = start of the final persistent snow
    cover of the winter (SNODP > ONSET_DEPTH_THRESHOLD): the date snow began
    accumulating and then stayed (apart from gaps shorter than ONSET_GAP_DAYS)
    until spring melt-out. If the pack melts out completely in November and
    re-accumulates in December until June, the December date is returned.
    The season searched runs from ONSET_SEARCH_START of year Y to
    ONSET_SEARCH_END of year Y+1.
    The CSV gives the ensemble-mean onset date for every year."""
    years = np.arange(FIRST_YEAR, LAST_YEAR + 1)
    onset = np.full((len(nc_files), len(years)), np.datetime64('NaT'), dtype='datetime64[ns]')
    run_names = []
    provisional_years = set()

    t0 = time.time()
    for r, fp in enumerate(nc_files):
        run_names.append(os.path.basename(fp).replace('.nc', ''))
        try:
            with xr.open_dataset(fp) as ds:
                s = ds['SNODP'].load()
        except Exception as e:
            print(f"  Error reading {fp}: {e}")
            continue

        t = pd.DatetimeIndex(s.time.values)
        depth = pd.Series(s.values.reshape(len(t)), index=t).sort_index()
        dt_s = np.median(np.diff(depth.index.values)) / np.timedelta64(1, 's')
        steps_per_day = 86400 / dt_s
        gap_n = max(1, int(round(ONSET_GAP_DAYS * steps_per_day)))          # steps
        min_n = max(1, int(round(ONSET_MIN_SEASON_DAYS * steps_per_day)))   # steps

        for j, y in enumerate(years):
            w = depth[f'{y}-{ONSET_SEARCH_START}':f'{y + 1}-{ONSET_SEARCH_END} 23:59']
            # skip years the run does not cover (needs data up to at least 31 Oct)
            if len(w) < gap_n or w.index[-1] < pd.Timestamp(f'{y}-10-31'):
                continue
            # run ends before the end of the winter -> onset may still change
            if w.index[-1] < pd.Timestamp(f'{y + 1}-{ONSET_PROVISIONAL_BEFORE}'):
                provisional_years.add(y)
            snow = (w > ONSET_DEPTH_THRESHOLD).values
            start = _final_snow_period_start(snow, w.index, gap_n, min_n)
            if start is not None:
                onset[r, j] = start.to_datetime64()

        if (r + 1) % max(1, len(nc_files) // 10) == 0 or r == len(nc_files) - 1:
            print(f"  Onset: processed {r + 1}/{len(nc_files)} runs ({time.time() - t0:.0f}s)")

    if provisional_years:
        print(f"\nWARNING: simulation ends before {ONSET_PROVISIONAL_BEFORE} of the following "
              f"year for {sorted(provisional_years)}; those onset dates are provisional "
              f"(they assume the snow present at the end of the run persists).")

    # Summary across members
    rows = []
    for j, y in enumerate(years):
        col = onset[:, j]
        valid = col[~np.isnat(col)]
        if valid.size == 0:
            rows.append((y, 0, '', np.nan, np.nan, '', ''))
            continue
        ts = pd.DatetimeIndex(valid)
        jan1 = pd.Timestamp(f'{y}-01-01')
        days = np.asarray((ts - jan1) / pd.Timedelta('1D'))        # fractional days since 1 Jan of year Y
        mean_days = days.mean()
        mean_date = jan1 + pd.Timedelta(days=mean_days)
        rows.append((y, valid.size, mean_date.strftime('%Y-%m-%d'),
                     mean_days + 1.0,                              # mean day of year (can exceed 365)
                     days.std(),
                     ts.min().strftime('%Y-%m-%d'), ts.max().strftime('%Y-%m-%d')))

    table = pd.DataFrame(rows, columns=['year', 'n_members_with_onset', 'mean_onset_date',
                                        'mean_onset_doy', 'std_days',
                                        'earliest_onset', 'latest_onset'])
    table['mean_onset_doy'] = table['mean_onset_doy'].round(1)
    table['std_days'] = table['std_days'].round(1)
    table.to_csv(csv_path, index=False)
    print(f"\nSnow onset (ensemble mean; start of final persistent snow cover, "
          f"depth > {ONSET_DEPTH_THRESHOLD} m, gaps < {ONSET_GAP_DAYS} d ignored, "
          f"min. duration {ONSET_MIN_SEASON_DAYS} d)")
    print(table.to_string(index=False))
    print(f"Saved {csv_path}")

    if members_csv_path:
        members = pd.DataFrame(onset, index=run_names, columns=years)
        members = members.apply(lambda c: c.dt.strftime('%Y-%m-%d'))
        members.index.name = 'simulation'
        members.to_csv(members_csv_path)
        print(f"Saved {members_csv_path}")
    return table


# ------------------------------------------------------------------------- #
# Plotting helpers
# ------------------------------------------------------------------------- #
def make_cmap_norm(n_runs):
    cmap = plt.cm.CMRmap_r
    colors = [cmap(i) for i in range(cmap.N)]
    colors[0] = (1.0, 1.0, 1.0, 1.0)        # 0 % -> white, as in Figure 7
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colors, cmap.N)
    norm = mpl.colors.BoundaryNorm(np.linspace(0, 100, n_runs + 1), cmap.N)
    return cmap, norm


def teardrop_marker():
    t = np.linspace(0, 2 * np.pi, 60)
    verts = np.column_stack([np.sin(t) * np.sin(t / 2), np.cos(t)])
    verts[:, 1] = verts[:, 1] * 0.9
    return Path(verts, closed=True)


def add_event_markers(ax, years):
    """Wet-snow detections (drops) and die-off years (stars) stacked from the
    bottom of each column; positions are in axes fraction so they work for
    both the normalized and the metre-based figure."""
    trans = ax.get_xaxis_transform()
    drop = teardrop_marker()
    step, y0 = 0.045, 0.03
    for y in years:
        k = WET_SNOW_DETECTIONS.get(y, 0)
        for i in range(k):
            ax.plot(y, y0 + i * step, marker=drop, ms=7, color='k', mec='none',
                    ls='none', transform=trans, zorder=6, clip_on=False)
        if y in ROS_DIEOFF_YEARS:
            ax.plot(y, y0 + k * step, marker='*', ms=9, color='k', mec='none',
                    ls='none', transform=trans, zorder=6, clip_on=False)
    handles = [Line2D([], [], marker=drop, ms=7, color='k', ls='none', label='Wet snow detection'),
               Line2D([], [], marker='*', ms=9, color='k', ls='none',
                      label='ROS event associated with ungulate die-off')]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.0),
              ncol=2, frameon=False)


def plot_ensemble(ds, redistribute, outfile):
    years = ds.year.values
    prob = np.ma.masked_invalid(ds['prob_high_density'].values)   # (year, bin)
    N = prob.shape[1]
    frac_edges = np.linspace(0, 1, N + 1)
    n_runs = ds.sizes['model_run']

    mean_depth = np.nanmean(ds['total_snow_depth'].values, axis=0)  # (year,)

    # Drop leading/trailing years with no data at all (e.g. a simulation that
    # ends before 31 Oct of the last year) so the left/right margins are even.
    has_data = ~np.ma.getmaskarray(prob).all(axis=1)
    idx = np.where(has_data)[0]
    sl = slice(idx[0], idx[-1] + 1)
    years, prob, mean_depth = years[sl], prob[sl], mean_depth[sl]

    plt.style.use('ggplot')
    fig, ax = plt.subplots(figsize=(11, 6))
    cmap, norm = make_cmap_norm(n_runs)

    half = COLUMN_WIDTH / 2
    for j, yr in enumerate(years):
        if prob[j].mask.all():
            continue
        if redistribute:
            if not np.isfinite(mean_depth[j]) or mean_depth[j] <= 0:
                continue
            y_edges = frac_edges * mean_depth[j]
        else:
            y_edges = frac_edges
        ax.pcolormesh([yr - half, yr + half], y_edges, prob[j][:, None],
                      cmap=cmap, norm=norm, shading='flat', zorder=2)
        # Steelblue outline around the whole column (same colour as the mean-depth
        # line in Figure 7). Columns touch, so neighbours share one edge.
        ax.add_patch(Rectangle((yr - half, 0), COLUMN_WIDTH, y_edges[-1], fill=False,
                               edgecolor=OUTLINE_COLOR, linewidth=0.8, zorder=4))

    if redistribute:
        # Background silhouette + dashed line = ensemble-mean snow depth
        ok = np.isfinite(mean_depth)
        ax.set_ylim(0, np.nanmax(mean_depth) * 1.05)
        ax.set_ylabel('Ensemble-mean snow depth on 31 October (m)', fontsize=12)
        for yr, d in zip(years[ok], mean_depth[ok]):
            ax.fill_between([yr - half, yr + half], [d, d], 0, color='whitesmoke', zorder=1.5)
    else:
        ax.set_ylim(0, 1)
        ax.set_ylabel('Normalized snow depth', fontsize=12)

    # Margins: SIDE_MARGIN bar-widths of empty space on both sides.
    ax.set_xlim(years[0] - COLUMN_WIDTH / 2 - SIDE_MARGIN * COLUMN_WIDTH,
                years[-1] + COLUMN_WIDTH / 2 + SIDE_MARGIN * COLUMN_WIDTH)
    # Major ticks (labelled) every TICK_STEP_YEARS bars starting at the first
    # year shown (1985, 1988, ...); minor ticks (unlabelled) under every bar.
    major = years[::TICK_STEP_YEARS]
    ax.set_xticks(major)
    ax.set_xticklabels([str(y) for y in major])
    ax.set_xticks(years, minor=True)
    ax.tick_params(axis='x', which='major', length=7, width=1, color='black')
    ax.tick_params(axis='x', which='minor', length=3, width=0.8, color='black')
    ax.set_xlabel('Year', fontsize=12)
    ax.grid(True, which='major', zorder=1)   # grid follows major ticks only
    ax.grid(False, which='minor')
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color('black')

    add_event_markers(ax, years)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.04, ticks=np.linspace(0, 100, 11),
                      boundaries=norm.boundaries, spacing='proportional', format='%1i')
    cb.set_label(f'Probability of high density snow (> {DENSITY_THRESHOLD} kg m$^{{-3}}$) (%)', size=12)
    for s in cb.ax.spines.values():
        s.set_visible(True)
        s.set_color('black')

    fig.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved {outfile}")
    return fig


# ------------------------------------------------------------------------- #
# Main
# ------------------------------------------------------------------------- #
if __name__ == '__main__':
    nc_files = None
    if os.path.exists(data_filename) and not args.recompute:
        print(f"Using existing {data_filename} (use --recompute to rebuild)")
        ds = xr.open_dataset(data_filename)
    else:
        nc_files, base = find_member_files()
        print(f"Found {len(nc_files)} members")
        ds = build_dataset(nc_files, base)

    # Quick table: probability that the basal layer is a locked pasture
    basal = ds['prob_high_density'].isel(depth_fraction=0).values
    print("\nYear  members_with_snow  P(basal layer > %d kg/m3) [%%]" % DENSITY_THRESHOLD)
    for y, n, p in zip(ds.year.values, ds['n_members_with_snow'].values, basal):
        print(f"{y}  {int(n):>17d}  {p:6.1f}" if np.isfinite(p) else f"{y}  {int(n):>17d}     n/a")

    # Bottom 20 mm values, as drawn in the redistributed figure
    print_basal_values(ds, thickness_mm=BASAL_THICKNESS_MM,
                       years_of_interest=BASAL_YEARS_OF_INTEREST,
                       csv_path=os.path.join(result_folder, 'basal_20mm_probabilities.csv'))

    # Ensemble-mean snow onset date for every year (re-reads SNODP from all members)
    if not args.no_onset:
        if nc_files is None:
            nc_files, _ = find_member_files()
        compute_snow_onset_table(
            nc_files,
            csv_path=os.path.join(result_folder, 'snow_onset_ensemble_mean.csv'),
            members_csv_path=os.path.join(result_folder, 'snow_onset_all_members.csv'))

    plot_ensemble(ds, redistribute=True,
                  outfile=os.path.join(result_folder, 'locked_pastures_oct31_redistributed.png'))
    plt.show()
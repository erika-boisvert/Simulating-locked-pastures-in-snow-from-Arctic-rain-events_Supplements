import datetime
import glob
import os
import pdb
import shutil
import subprocess
import sys

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import collections
from matplotlib.collections import PatchCollection
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.colors import Normalize
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# --- PARAMETERS & ARGUMENTS ---
if len(sys.argv) != 3:
    raise ValueError('Requires one argument: <name of the experiment>')

sim_name = sys.argv[-2]
event_year = sys.argv[-1]

os.chdir(os.path.join(r'/home/ebv000/experiments', sim_name, 'output'))    # Navigate to the output folder
nc_path = r'./out_svs2.nc'    #output netcdf standard name
# nc_path = r'./all_outputs/CaSR32R_norain_HOLDB92_COMPS14_FALLR21_DRIFTR21F.nc'
ds = xr.open_dataset(nc_path)

date = ds.time[-1].values
year = date.astype('datetime64[Y]').astype(int) + 1970

date_beg = pd.to_datetime(str(event_year) + '-07-02')
# date_beg = ds.time.values.min() + np.timedelta64(1, 'D')
date_end = pd.to_datetime(str(int(event_year)+1) + '-06-30')

mod = ds.loc[dict(time=slice(date_beg,date_end))]

output_name = sim_name + '_' + str(event_year) + '_plot.png'

# --- FUNCTIONS ---
def import_bulk_sim(mod, exp='opn'):
    # exp is either "opn" or "for" depending on if we want to analyse the snowpack in the open or in a forest
    if exp == 'opn':
        swe = mod['SNOMA'].to_dataframe('swe')
        sd = mod['SNODP'].to_dataframe('sd')
        Ts = mod['TSNO_SURF'].to_dataframe('Ts')
        try:
            snow_can = mod['SNOMA'].to_dataframe('snow_can')
            sim = pd.concat([swe, sd, Ts, snow_can ], axis=1)
        except:
            sim = pd.concat([swe, sd, Ts ], axis=1)
        sim['den'] = sim.swe.values / sim.sd.values
        Thickness = mod['SNOMA_ML'].values / mod['SNODEN_ML'].values
        T_snow = mod['TSNOW_ML'].values

    elif exp == 'for':
        swe = mod['SNVMA'].to_dataframe('swe')
        sd = mod['SNVDP'].to_dataframe('sd')
        Ts = mod['TSNV_SURF'].to_dataframe('Ts')
        snow_can = mod['SNCMA'].to_dataframe('snow_can')
        sim = pd.concat([swe, sd, Ts, snow_can ], axis=1)
        sim['den'] = sim.swe.values / sim.sd.values
        Thickness = mod['SNOMAV_ML'].values / mod['SNODENV_ML'].values
        T_snow = mod['TSNOWV_ML'].values

    Tmean = []
    for i in range(len(T_snow)):

        T = T_snow[i][Thickness[i] > 0]
        thick = Thickness[i][Thickness[i] > 0]
        if len(T) > 0:
            Tmean.append(np.average(T, weights = thick))
        else:
            Tmean.append(np.nan)
    Tmean = np.array(Tmean)

    ci = 2.1e-3 # MJ/kg/C
    CC = (273.15-Tmean) * ci * sim['den'].values * sim['sd'].values
    sim['CC'] = CC

    return sim

def import_snow_layers_sim(mod, exp='opn'):
    # exp is either "opn" or "for" depending on if we want to analyse the snowpack in the open or in a forest
    # if exp == 'opn':
    mod['Thickness'] = mod.variables['SNOMA_ML']/mod.variables['SNODEN_ML']
    df = mod[['SNODEN_ML','SNOMA_ML','SNODP','SNODOPT_ML','TSNOW_ML','SNOTYPE_ML', 'Thickness', 'WSNOW_ML', 'SNOAGE_ML']].to_dataframe()
    df['SSA'] = df['SNODOPT_ML'].apply(lambda x : 6./(x*917) if x > 0 else np.nan)
    # elif exp == 'for':
    #     df = mod[['SNODENV_ML','SNOMAV_ML','SNVDP','SNODOPTV_ML','TSNOWV_ML','SNOTYPEV_ML', 'WSNOW_ML']].to_dataframe()
    #     df['Thickness'] = df[['SNODENV_ML','SNOMAV_ML']].apply(lambda x : x.iloc[1] / x.iloc[0], axis = 1)
    #     df['SSA'] = df['SNODOPTV_ML'].apply(lambda x : 6./(x*917) if x > 0 else np.nan)
    df.columns = ['rho_ML','SWE_ML','Depth','Dgrain_ML','T_ML','Gtype_ML','Thickness_ML', 'WSNOW_ML', 'SNOAGE_ML', 'SSA_ML']
    return df

def import_soil_layers_sim(mod):
    # exp is either "opn" or "for" depending on if we want to analyse the snowpack in the open or in a forest
    df = mod[['WSOIL','ISOIL','TPSOIL']].to_dataframe()
    df.columns = ['Water_content','Ice_content','Temperature_K']
    return df

def get_season_profile(df, variable, nlayers=20):
    time = df.index.get_level_values('time').unique()
    # print(time)
    # assert False

    dz = df['Thickness_ML'].values
    dz = np.reshape(dz, (-1, nlayers))
    dz[np.isnan(dz)] = 0

    value = df[variable].values
    value = np.reshape(value, (-1,nlayers))
    value[np.isnan(value)] = 0

    top_y = np.cumsum(dz[:, ::-1], axis=1)[:, ::-1].ravel()
    bottom_y = top_y - dz.ravel()

    left_x = np.ones(shape=dz.shape, dtype='int')
    left_x = np.cumsum(left_x, axis=0).ravel()
    # left_x = np.arange(dz.shape[0])
    right_x = left_x + 9
    # right_x = left_x - 1
    mid_x = (left_x + right_x)/2
    mid_x = mid_x[::nlayers]

    vertices = np.zeros(shape=(right_x.shape[0], 4, 2))
    vertices[:, 0, 0] = right_x
    vertices[:, 0, 1] = bottom_y
    vertices[:, 1, 0] = right_x
    vertices[:, 1, 1] = top_y
    vertices[:, 2, 0] = left_x
    vertices[:, 2, 1] = top_y
    vertices[:, 3, 0] = left_x
    vertices[:, 3, 1] = bottom_y
    vertices = vertices[(dz > 0).ravel()]
    return time, vertices, value, dz, mid_x

def plot_2Dprofile(ax, times, df_sim, variable):

    n = len(times) 
    if variable == 'Gtype_ML':
        cmap = grain_colormap
        bounds = np.linspace(-0.5, 14.5, 16)
        norm = BoundaryNorm(bounds, cmap.N)
    else :
        cmap = 'jet'

    max_height = 0
    for j, time in enumerate(times):
        Type = df_sim.loc[time][variable].values
        Type[np.isnan(Type)] = np.nan
        
        height_max = 0
        hgt, hgt_norm, value = pit_sim(df_sim, time, variable)
        
        hgt_bot = hgt[1::2]
        hgt_top = hgt[::2]
        max_height = max([max_height,max(hgt_top)])
    
        # Loop over data points; create box from snowtype  at each layer
        layers_sim = [Rectangle((j/n, h_bot), (j+1)/n - j/n, h_top-h_bot,edgecolor=None) for h_bot, h_top in zip(list(hgt_bot),list(hgt_top))]
        height_max = np.max([height_max, np.max(hgt_top)])
    
        # Create patch collection with specified colour/alpha
        if variable == 'Gtype_ML':
            pc = PatchCollection(layers_sim, array = Type, cmap=cmap, norm=norm)#, edgecolor='black')
        else :
            pc = PatchCollection(layers_sim, array = Type, cmap=cmap)#, edgecolor='black')

        # Add collection to axes
        ax.add_collection(pc)
        plt.axvline(j/n, color = 'grey')
        
    ax.set_xticks(np.linspace(0.5/n, 1.-0.5/n, n), times, rotation = 45)
    ax.set_xlim(0,1)
    ax.set_ylim(0,max_height + 0.1)
    cbar = plt.colorbar(pc, ax=ax)

    if variable == 'Gtype_ML':
        labels = MEPRA_labels
        cbar.set_ticks(np.arange(np.shape(labels)[0]))
        cbar.ax.set_yticklabels(labels)
    else :
        cbar.ax.set_title(variable)
    return (pc, height_max)

def get_plot_coord(plot_no):
    if plot_no == 0:
        i,j = 0,0
    elif plot_no == 1:
        i,j = 1,0
    elif plot_no == 2:
        i,j = 0,1
    elif plot_no == 3:
        i,j = 1,1
    return i, j

def format_ticks(x, pos):
    x = int(x)
    x = max(min(len(time) - 1, x), 0)
    if hasattr(time[0], 'strftime'):
        return time[x].strftime('%Y-%m-%d')
    else:
        return time[x]

def weighted_average(df, density, thickness):
    return sum(df[thickness] * df[density]) / df[thickness].sum()

def basal_dense_weighted_density(df_single_time, density_threshold=350,
                                   density_col='rho_ML', thickness_col='Thickness_ML'):
    """Given the layers for ONE timestep (indexed by snow_layer, ascending --
    index 0 = surface/top, highest valid index = ground contact, matching the
    SNOMA_ML/SNODEN_ML convention used elsewhere), find the contiguous set of
    layers starting at the ground and moving upward that all exceed
    `density_threshold`, and return the thickness-weighted average density of
    just that basal region.

    Returns (basal_density, basal_thickness); both np.nan/0.0 if the
    ground-contact layer itself isn't dense enough (no basal layer at all).
    """
    d = df_single_time.sort_index()  # ensure ascending snow_layer order (top -> ground)
    thickness_vals = d[thickness_col].values
    density_vals = d[density_col].values

    valid = thickness_vals > 0
    if not valid.any():
        return np.nan, 0.0

    ground_idx = np.where(valid)[0].max()  # deepest layer with actual snow

    dense_layers = []
    for idx in range(ground_idx, -1, -1):  # walk upward starting at the ground
        if thickness_vals[idx] <= 0 or density_vals[idx] <= density_threshold:
            break
        dense_layers.append(idx)

    if not dense_layers:
        return np.nan, 0.0

    basal_thickness = thickness_vals[dense_layers].sum()
    basal_density = np.average(density_vals[dense_layers], weights=thickness_vals[dense_layers])
    return basal_density, basal_thickness

# --- DATA DICTIONARIES ---
## From proplotter (Snowtool)
# The list of SURFEX/MEPRA grain shapes associated to its common names

MEPRA_dict = {'PP': 0,      # fr
              'PP+DF': 1,   # fr_lb
              'DF': 2,      # lb
              'DF+RG': 3,   # lb_fin
              'DF+FC': 4,   # lb_ang
              'PPgp': 5,    # roul
              'RG': 6,      # fin
              'MF+RG': 7,   # fin_ar
              'RG+FC': 8,   # fin_ang
              'FC': 9,      # pl
              'FC+DH': 10,  # pl_gob
              'DH': 11,     # gob
              'MF': 12,     # gel
              'MF+DH': 13,  # gob_fon
              'MF+FC': 14   # ron_ang
              }

MEPRA_labels = ['PP', 'PP+DF', 'DF', 'DF+RG', 'DF+FC', 'PPgp', 'RG', 'MF+RG', 'RG+FC', 'FC', 'FC+DH', 'DH', 'MF', 'MF+DH', 'MF+FC']

coloring = {'PP': np.array([0, 255, 0]) / 255.0,
            'MM': np.array([255, 215, 0]) / 255.0,
            'DF': np.array([34, 139, 34]) / 255.0,
            'RG': np.array([255, 182, 193]) / 255.0,
            'FC': np.array([173, 216, 230]) / 255.0,
            'DH': np.array([0, 0, 255]) / 255.0,
            'SH': np.array([250, 0, 255]) / 255.0,
            'MF': np.array([255, 0, 0]) / 255.0,
            'MFcr': np.array([255, 255, 255]) / 255.0,
            'IF': np.array([0, 255, 255]) / 255.0,
            'NO': np.array([200, 200, 200]) / 255.0}

color_grain = {MEPRA_dict['PP']: coloring['PP'],
               MEPRA_dict['PP+DF']: 0.5 * (coloring['PP'] + coloring['DF']),
               MEPRA_dict['DF']: coloring['DF'],
               MEPRA_dict['DF+RG']: 0.5 * (coloring['DF'] + coloring['RG']),
               MEPRA_dict['DF+FC']: 0.5 * (coloring['DF'] + coloring['FC']),
               MEPRA_dict['PPgp']: [0, 0, 0],
               MEPRA_dict['RG']: coloring['RG'],
               MEPRA_dict['MF+RG']: 0.5 * (coloring['MF'] + coloring['RG']),
               MEPRA_dict['RG+FC']: 0.5 * (coloring['RG'] + coloring['FC']),
               MEPRA_dict['FC']: coloring['FC'],
               MEPRA_dict['FC+DH']: 0.5 * (coloring['FC'] + coloring['DH']),
               MEPRA_dict['DH']: coloring['DH'],
               MEPRA_dict['MF']: coloring['MF'],
               MEPRA_dict['MF+DH']: 0.5 * (coloring['MF'] + coloring['DH']),
               MEPRA_dict['MF+FC']: 0.5 * (coloring['MF'] + coloring['FC'])}

grain_colormap = LinearSegmentedColormap.from_list("custom", [[i / 14., color_grain[i]] for i in range(15)])


# --- PREPARE THE DATA ---
df_bulk = import_bulk_sim(mod)
df_snowlayer = import_snow_layers_sim(mod)
df_soillayer = import_soil_layers_sim(mod)

# --- INSET SETTINGS (matched to normalized_ensemble_350) ---
INSET_WIDTH, INSET_HEIGHT = 0.30, 0.65
zoom_x1, zoom_x2 = 2145, 2375  # same index-based zoom window as before

n_steps = len(mod.time)
step_idx = np.arange(1, n_steps + 1)  # matches mid_x's indexing convention
zoom_mask_snodp = (step_idx >= zoom_x1) & (step_idx <= zoom_x2)

snodp_vals = mod['SNODP'].values
zoom_ylim_max = np.nanmax(snodp_vals[zoom_mask_snodp]) if zoom_mask_snodp.any() else np.nan
zoom_ylim_max = zoom_ylim_max * 1.1 if np.isfinite(zoom_ylim_max) and zoom_ylim_max > 0 else 0.055


# --- GET THE DENSITY OF THE DENSE BASAL LAYER (density > threshold, contiguous from ground) ---
density_threshold = 350
sub = df_snowlayer[['Thickness_ML', 'rho_ML']][str(event_year) + '-10-31': str(event_year) +'-11-01']
first_time = sub.index.get_level_values('time').unique()[0]  # same single-snapshot behavior as before
snapshot = sub.xs(first_time, level='time')

basal_density_val, basal_thickness_val = basal_dense_weighted_density(
    snapshot, density_threshold=density_threshold
)
basal_density = int(basal_density_val) if not np.isnan(basal_density_val) else None


# --- PLOT SETUP ---
plt.style.use('ggplot')

fig = plt.figure(figsize=(18, 9))
gs = fig.add_gridspec(nrows = 2, ncols = 4, width_ratios=[1, 0.03, 1, 0.03], wspace = 0.25)

ax = np.empty((2, 2), dtype=object)
ax[0, 0] = fig.add_subplot(gs[0, 0])
ax[0, 1] = fig.add_subplot(gs[0, 2])
ax[1, 0] = fig.add_subplot(gs[1, 0])
ax[1, 1] = fig.add_subplot(gs[1, 2]) 

cax = {}
cax[(0, 0)] = fig.add_subplot(gs[0, 1]) 
cax[(0, 1)] = fig.add_subplot(gs[0, 3]) 
cax[(1, 0)] = fig.add_subplot(gs[1, 1]) 
cax[(1, 1)] = fig.add_subplot(gs[1, 3]) 


# fig, ax = plt.subplots(2, 2, figsize=(18, 9), sharex = 'all')
# fig.suptitle('Crocus SVS2 Output for %s from %s to %s'%(sim_name, np.str_(date_beg)[:10], np.str_(date_end)[:10]), fontsize = 15)

plot_no = 0

# --- COLOR SCHEMES ---
c_soil = ['goldenrod', 'darkkhaki', 'yellowgreen', 'green', 'darkcyan', 'purple', 'deeppink', 'darkorange', 'firebrick', 'chocolate', 'slategray', 'k'] #soil layers plot
cmap_dens = plt.cm.PuBu    # density plot
cmap_temp = plt.cm.coolwarm    # snowpack temperature plot
cmap_gtype = grain_colormap
cmap_water = plt.cm.BuPu


# --- PLOT DENSITY ---
variable = 'rho_ML'
time, vertices, values, dz, mid_x = get_season_profile(df_snowlayer, variable, nlayers = mod.snow_layer.shape[0]) # If the "time" output has another name, it should be changed in the fonction format_ticks()

#Wet snow detection period
# start1_i = np.where(time == np.datetime64('2003-10-02T12:00:00'))[0][0]
# end1_i = np.where(time == np.datetime64('2003-10-03T00:00:00'))[0][0]
# start2_i = np.where(time == np.datetime64('2003-10-04T00:00:00'))[0][0]
# end2_i = np.where(time == np.datetime64('2003-10-05T12:00:00'))[0][0]


i,j = get_plot_coord(plot_no)
# ax[i,j].plot(df_snowlayer.xs(0, level = 'snow_layer')['Depth'].index, df_snowlayer.xs(0, level = 'snow_layer')['Depth'], linewidth=1.0, c = 'k', alpha = 0.5)
ax[i,j].plot(mid_x, mod['SNODP'].values,linewidth=1.0, c = 'k', alpha = 0.5, zorder = 4)

norm = Normalize(vmin=150, vmax=550)
rect = collections.PolyCollection(
    vertices, 
    array=values[(dz > 0)].ravel(),
    cmap=cmap_dens, 
    edgecolors='none',
    norm = norm, 
    zorder = 3
)
ax[i,j].add_collection(rect)
ax[i,j].autoscale_view()

cb1 = fig.colorbar(rect, cax=cax[(i,j)], orientation='vertical')
cb1.ax.set_title('ρ (kg m$^{-3}$)', fontsize=10)

# cbar = plt.colorbar(rect, ax=ax[i,j])
# cbar.ax.set_title('Density (kg/m³)', fontsize = 8)
# fig.colorbar(cbar, ax=ax[i,j], orientation='vertical')

formatter = ticker.FuncFormatter(format_ticks)
ax[i,j].xaxis.set_major_formatter(formatter)
plt.setp(ax[i,j].xaxis.get_majorticklabels(), size='small')
plt.xticks(rotation = 70)


ax[i,j].set_title('a) Snowpack Density Evolution', loc = 'left')
ax[i,j].set_ylabel('Snow depth (m)', fontsize = 12)

if basal_density is not None:
    ax[i,j].text(4500, 0.010, str(basal_density) + ' kg m$^{-3}$', fontsize = 14, weight = 'semibold', zorder = 6)
else:
    ax[i,j].text(5000, 0.026, 'No dense basal layer', fontsize = 12, weight = 'bold', zorder = 6)

pos = ax[i, j].get_position()
ax[i, j].set_position([
    pos.x0 + 0.02,  # shift right
    pos.y0,
    pos.width,
    pos.height
])

ax[i,j].set_ylim([0, 0.5])

## ADD INSET
inset_ax = ax[i,j].inset_axes(
   [0.0, 1.0 - INSET_HEIGHT, INSET_WIDTH, INSET_HEIGHT],  # matches normalized_ensemble_350
    xlim=[zoom_x1, zoom_x2], ylim=[0.0, zoom_ylim_max],
    xticklabels=[], yticklabels=[]
)

for axi in ax[i,j], inset_ax:
    rect1 = collections.PolyCollection(vertices, array=values[(dz > 0)].ravel(), 
                                      cmap=cmap_dens, edgecolors='none', norm = norm, zorder = 3)
    axi.add_collection(rect1)  
    axi.plot(mid_x, mod['SNODP'].values,linewidth=1.0, c = 'k', alpha = 0.5, zorder = 4) 
    # axi.axvspan(mid_x[start1_i], mid_x[end1_i], color='red', alpha=0.16, zorder = 1)
    # axi.axvspan(mid_x[start2_i], mid_x[end2_i], color='red', alpha=0.16, zorder = 1)


ax[i,j].indicate_inset_zoom(inset_ax, edgecolor="k", linewidth=1.0)
for spine in inset_ax.spines.values():
    spine.set_edgecolor("k")       # frame color 
plot_no += 1


# --- WATER CONTENT ---
variable = 'WSNOW_ML'
time, vertices, values, dz, mid_x = get_season_profile(df_snowlayer, variable, nlayers = mod.snow_layer.shape[0]) # If the "time" output has another name, it should be changed in the fonction format_ticks()

#Wet snow detection period
# start1_i = np.where(time == np.datetime64('2003-10-02T12:00:00'))[0][0]
# end1_i = np.where(time == np.datetime64('2003-10-03T00:00:00'))[0][0]
# start2_i = np.where(time == np.datetime64('2003-10-04T00:00:00'))[0][0]
# end2_i = np.where(time == np.datetime64('2003-10-05T12:00:00'))[0][0]

i,j = get_plot_coord(plot_no)
# ax[i,j].plot(df_snowlayer.index.get_level_values("time").unique(), df_snowlayer.xs(0, level = 'snow_layer')['Depth'], linewidth=1.0, c = 'k', alpha = 0.5)
ax[i,j].plot(mid_x, mod['SNODP'].values,linewidth=1.0, c = 'k', alpha = 0.5, zorder = 4)

norm = Normalize(vmin=0, vmax=0.0017)
rect = collections.PolyCollection(
    vertices, 
    array=values[(dz > 0)].ravel(),
    cmap=cmap_water, 
    edgecolors='none',
    zorder = 3,
    norm = norm
)
ax[i,j].add_collection(rect)
ax[i,j].autoscale_view()

cb2 = fig.colorbar(rect, cax=cax[(i,j)], orientation='vertical')
cb2.ax.set_title('LWC (kg m$^{-2}$)', fontsize=10)

# cbar = plt.colorbar(rect, ax=ax[i,j])
# cbar.ax.set_title('LWC (kg/m²)',  fontsize = 8)

formatter = ticker.FuncFormatter(format_ticks)
ax[i,j].xaxis.set_major_formatter(formatter)
plt.setp(ax[i,j].xaxis.get_majorticklabels(), size='small')
plt.xticks(rotation = 90)

ax[i,j].set_title('c) Snowpack Liquid Water Content (LWC)', loc = 'left')
ax[i,j].set_ylabel('Snow depth (m)', fontsize = 12)

pos = ax[i, j].get_position()
ax[i, j].set_position([
    pos.x0 + 0.02,  # shift right
    pos.y0,
    pos.width,
    pos.height
])

## ADD INSET
ax[i,j].set_ylim([0, 0.5])

inset_ax = ax[i,j].inset_axes(
   [0.0, 1.0 - INSET_HEIGHT, INSET_WIDTH, INSET_HEIGHT],  # matches normalized_ensemble_350
    xlim=[zoom_x1, zoom_x2], ylim=[0.0, zoom_ylim_max],
    xticklabels=[], yticklabels=[]
)

for axi in ax[i,j], inset_ax:
    rect1 = collections.PolyCollection(vertices, array=values[(dz > 0)].ravel(),
                                       cmap=cmap_water, edgecolors='none', norm = norm, zorder = 3)
    axi.add_collection(rect1)  
    axi.plot(mid_x, mod['SNODP'].values,linewidth=1.0, c = 'k', alpha = 0.5,  zorder = 3) 
    # axi.axvspan(mid_x[start1_i], mid_x[end1_i], color='red', alpha=0.16, zorder = 1)
    # axi.axvspan(mid_x[start2_i], mid_x[end2_i], color='red', alpha=0.16, zorder = 1)

ax[i,j].indicate_inset_zoom(inset_ax, edgecolor="k", linewidth=1.0)
for spine in inset_ax.spines.values():
    spine.set_edgecolor("k")       # frame color

plot_no += 1

# --- PLOT GRAIN TYPE ---
variable = 'Gtype_ML'

time, vertices, Gtype, dz, mid_x = get_season_profile(df_snowlayer, variable, nlayers = mod.snow_layer.shape[0]) # If the "time" output has another name, it should be changed in the fonction format_ticks()

#Wet snow detection period
# start1_i = np.where(time == np.datetime64('2003-10-02T12:00:00'))[0][0]
# end1_i = np.where(time == np.datetime64('2003-10-03T00:00:00'))[0][0]
# start2_i = np.where(time == np.datetime64('2003-10-04T00:00:00'))[0][0]
# end2_i = np.where(time == np.datetime64('2003-10-05T12:00:00'))[0][0]


bounds = np.linspace(-0.5, 14.5, 16)
norm = BoundaryNorm(bounds, cmap_gtype.N)
vmin = -0.5
vmax = 14.5
extend = 'neither'

i,j = get_plot_coord(plot_no)
# ax[i,j].plot(df_snowlayer.index.get_level_values("time").unique(), df_snowlayer.xs(0, level = 'snow_layer')['Depth'], linewidth=1.0, c = 'k', alpha = 0.5)
ax[i,j].plot(mid_x, mod['SNODP'].values, linewidth=1.0, c = 'k', alpha = 0.5,  zorder = 4)

rect = collections.PolyCollection(vertices, array=Gtype[(dz > 0)].ravel(),
                                    cmap=cmap_gtype, norm=norm, edgecolors='none', zorder = 3)
rect.set_clim(vmin, vmax)
ax[i,j].add_collection(rect)
ax[i,j].autoscale_view()


cb3 = fig.colorbar(rect, cax=cax[(i,j)], orientation='vertical', extend = extend)
labels = MEPRA_labels
cb3.set_ticks(np.arange(np.shape(labels)[0]))
cb3.ax.set_yticklabels(labels)


formatter = ticker.FuncFormatter(format_ticks)
ax[i,j].xaxis.set_major_formatter(formatter)
plt.setp(ax[i,j].xaxis.get_majorticklabels(), size='small')
plt.xticks(rotation = 90)

ax[i,j].set_title('b) Grain Type', loc = 'left')
ax[i,j].set_ylabel('Snow depth (m)', fontsize = 12)

pos = ax[i, j].get_position()
ax[i, j].set_position([
    pos.x0 + 0.02,  # shift right
    pos.y0,
    pos.width,
    pos.height
])

## ADD INSET
ax[i,j].set_ylim([0, 0.5])

inset_ax = ax[i,j].inset_axes(
   [0.0, 1.0 - INSET_HEIGHT, INSET_WIDTH, INSET_HEIGHT],  # matches normalized_ensemble_350
    xlim=[zoom_x1, zoom_x2], ylim=[0.0, zoom_ylim_max],
    xticklabels=[], yticklabels=[]
)
for axi in ax[i,j], inset_ax:
    rect1 = collections.PolyCollection(vertices, array=Gtype[(dz > 0)].ravel(),
                                    cmap=cmap_gtype, norm=norm, edgecolors='none')
    axi.add_collection(rect1)  
    axi.plot(mid_x, mod['SNODP'].values,linewidth=1.0, c = 'k', alpha = 0.5,  zorder = 4) 
    # axi.axvspan(mid_x[start1_i], mid_x[end1_i], color='red', alpha=0.16, zorder = 1)
    # axi.axvspan(mid_x[start2_i], mid_x[end2_i], color='red', alpha=0.16, zorder = 1)

ax[i,j].indicate_inset_zoom(inset_ax, edgecolor="k", linewidth=1.0)
for spine in inset_ax.spines.values():
    spine.set_edgecolor("k")       # frame color


plot_no += 1


# --- SOIL LAYER TEMPERATURE ---
# variable = 'Temperature_K'

# i,j = get_plot_coord(plot_no)
# n_soillayer = mod.soil_layer.shape[0]

# for s_layer in range(n_soillayer):
#     soil_temp_values = df_soillayer.xs(s_layer, level = 'soil_layer')[variable]
#     if len(soil_temp_values.unique()) != 1:
#         ax[i,j].plot(df_soillayer.index.get_level_values("time").unique(), soil_temp_values, label=s_layer, c = c_soil[s_layer], linewidth=1)
#         # ax[i,j].plot(df_soillayer.index.get_level_values("time").unique(), soil_temp_values, label=s_layer, linewidth=1)
# ax[i,j].axhline(273.16, zorder = 1, color = 'k', linestyle = '--', alpha = 0.3)
# # ax[0,1].set_xlabel('Time')
# ax[i,j].set_ylabel('Temperature (K)')
# ax[i,j].set_title('Soil Layers Temperature')
# ax[i,j].legend().remove()
# ax[i,j].legend(bbox_to_anchor=(1.15, 0.5), loc = 'center right', ncol = 1, frameon = False)

# pos = ax[i, j].get_position()
# ax[i, j].set_position([
#     pos.x0 + 0.02,  # shift right
#     pos.y0,
#     pos.width,
#     pos.height
# ])

# fig.delaxes(cax[(i,j)])
# plot_no += 1

# --- SNOW LAYER TEMPERATURE ---
variable = 'T_ML'
cmap_temp = plt.cm.turbo

time, vertices, values, dz, mid_x = get_season_profile(df_snowlayer, variable, nlayers = mod.snow_layer.shape[0]) # If the "time" output has another name, it should be changed in the fonction format_ticks()

#Wet snow detection period
# start1_i = np.where(time == np.datetime64('2003-10-02T12:00:00'))[0][0]
# end1_i = np.where(time == np.datetime64('2003-10-03T00:00:00'))[0][0]
# start2_i = np.where(time == np.datetime64('2003-10-04T00:00:00'))[0][0]
# end2_i = np.where(time == np.datetime64('2003-10-05T12:00:00'))[0][0]


i,j = get_plot_coord(plot_no)
ax[i,j].plot(mid_x, mod['SNODP'].values,linewidth=1.0, c = 'k', alpha = 0.5, zorder = 4)

norm = Normalize(vmin=200, vmax=275)

rect = collections.PolyCollection(
    vertices, 
    array=values[(dz > 0)].ravel(),
    cmap=cmap_temp, 
    edgecolors='none',
    norm = norm,
    zorder = 3
)
ax[i,j].add_collection(rect)
ax[i,j].autoscale_view()

cb2 = fig.colorbar(rect, cax=cax[(i,j)], orientation='vertical')
cb2.ax.set_title('Temperature (K)', fontsize=10)

# cbar = plt.colorbar(rect, ax=ax[i,j])
# cbar.ax.set_title('LWC (kg/m²)',  fontsize = 8)

formatter = ticker.FuncFormatter(format_ticks)
ax[i,j].xaxis.set_major_formatter(formatter)
plt.setp(ax[i,j].xaxis.get_majorticklabels(), size='small')
plt.xticks(rotation = 90)

ax[i,j].set_title('d) Snowpack Temperature Profile', loc = 'left')
ax[i,j].set_ylabel('Snow depth (m)', fontsize = 12)

pos = ax[i, j].get_position()
ax[i, j].set_position([
    pos.x0 + 0.02,  # shift right
    pos.y0,
    pos.width,
    pos.height
])

## ADD INSET
ax[i,j].set_ylim([0, 0.5])

inset_ax = ax[i,j].inset_axes(
   [0.0, 1.0 - INSET_HEIGHT, INSET_WIDTH, INSET_HEIGHT],  # matches normalized_ensemble_350
    xlim=[zoom_x1, zoom_x2], ylim=[0.0, zoom_ylim_max],
    xticklabels=[], yticklabels=[]
)

for axi in ax[i,j], inset_ax:
    rect1 = collections.PolyCollection(vertices, array=values[(dz > 0)].ravel(),
                                       cmap=cmap_temp, edgecolors='none', zorder = 3)
    axi.add_collection(rect1)  
    axi.plot(mid_x, mod['SNODP'].values,linewidth=1.0, c = 'k', alpha = 0.5,  zorder = 4) 
    # axi.axvspan(mid_x[start1_i], mid_x[end1_i], color='red', alpha=0.16, zorder = 1)
    # axi.axvspan(mid_x[start2_i], mid_x[end2_i], color='red', alpha=0.16, zorder = 1)

ax[i,j].indicate_inset_zoom(inset_ax, edgecolor="k", linewidth=1.0)
for spine in inset_ax.spines.values():
    spine.set_edgecolor("k")       # frame color


# --- SAVE FIG ---
for key in cax:
    cax[key].tick_params(left=False, labelleft=False)
    pos = cax[key].get_position()
    cax[key].set_position([
        pos.x0 - 0.015,  # You can tweak this value (e.g. -0.015 for more shift)
        pos.y0,
        pos.width,
        pos.height
    ])

# fig.tight_layout()
plt.savefig(r'./' + output_name, bbox_inches='tight')

subprocess.call(['ln', '-sf', r'./' + output_name, r'/home/ebv000/public_html/output_figures'])

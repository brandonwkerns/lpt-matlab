import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import glob
import os

dir = '../../data/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/timeclusters'
out_dir = '../../plots/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/systems'
plt.close('all')


## Read in MJO LPTs
mjo_file = '/home/orca/bkerns/data/lpt/trmm3/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin2.txt'
clumps_file = '/home/orca/bkerns/data/lpt/trmm3/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'

D_mjo = np.loadtxt(mjo_file, skiprows=1)
mjo={}
mjo['year'] = D_mjo[:,0]
mjo['idx'] = D_mjo[:,1]
mjo['start'] = D_mjo[:,8]
mjo['end'] = D_mjo[:,9]

D_clumps = np.loadtxt(clumps_file, skiprows=1)
clumps={}
clumps['year'] = D_clumps[:,0]
clumps['idx'] = D_clumps[:,1]
clumps['clump_num'] = D_clumps[:,2]

os.system('mkdir -p ' + out_dir)


grand_total_count = 0
grand_total_lpt_count = 0


fout = (out_dir + '/lpt_mjo_system_tracks_20yrs__map_with_lpt.png')

## Set up figure
    
fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(211)
    
## Plot
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(llcrnrlon=0.,llcrnrlat=-50.,urcrnrlon=360.,urcrnrlat=50.,\
              resolution='l',projection='merc')
# draw coastlines, country boundaries, fill continents.
map.fillcontinents(color='wheat',zorder=1)
map.drawcoastlines(color='k',linewidth=0.7,zorder=10)
map.drawcountries(color='k',linewidth=0.5,zorder=10)
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary()
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30),zorder=15)
map.drawparallels(np.arange(-90,90,30),zorder=15)


for fn in sorted(glob.glob(dir+'/TIMECLUSTERS_*.rejoin2.nc')):

    # Skip RT 2018 for now.
    if "2018112721" in fn:
        continue
    
    print(fn, flush=True)
    year = int(fn[-32:-28])

    out_str = fn[-32:-11]

    clumps_this_year = {}
    clumps_this_year['year'] = clumps['year'][clumps['year'] == year]
    clumps_this_year['idx'] = clumps['idx'][clumps['year'] == year]
    clumps_this_year['clump_num'] = clumps['clump_num'][clumps['year'] == year]


                           

    ## Read
    DS = Dataset(fn)

    lon = DS['centroid_lon'][:]
    lat = DS['centroid_lat'][:]
    lptid = DS['id'][:]    
    time = DS['time'][:]
    area = DS['area'][:] / (1e5)

    mjo_count = 0
    lpt_count = 0
    
    for cnum in sorted(np.unique(clumps_this_year['clump_num'])):

        lpt_count += 1
        grand_total_lpt_count += 1

        ## mask out non-LPT systems.
        for idx in np.unique(clumps_this_year['idx'][clumps_this_year['clump_num'] == cnum]):

            grand_total_lpt_count += 1

            lon_this = lon[lptid == idx]
            lat_this = lat[lptid == idx]

            ## grey track in background.
            x, y = map(lon_this, lat_this)            
            h00 = map.plot(x,y,zorder=2, color=[0.5, 0.5, 0.5], linewidth=0.7)

            
            if np.logical_and(mjo['year'] == year, mjo['idx'] == idx).any():

                mjo_count += 1
                
                idx11 = int(mjo['start'][np.logical_and(mjo['year'] == year, mjo['idx'] == idx)][0])
                idx22 = int(mjo['end'][np.logical_and(mjo['year'] == year, mjo['idx'] == idx)][0])
                
                x, y = map(lon_this[idx11:idx22], lat_this[idx11:idx22])
                
                h0 = map.plot(x, y, color='yellowgreen', zorder=5, linewidth=0.7)
                h1 = map.plot(x[0], y[0], 'ko', markersize=2, zorder=1000)
                h2 = map.plot(x[-1], y[-1], 'rx', markersize=5, zorder=1000, markeredgewidth=1.5)
            
    print(mjo_count)
    grand_total_count += mjo_count


leg = plt.legend((h00[0], h0[0], h1[0], h2[0],),('LPT','MJO LPT', 'MJO Start','MJO End',), loc=(0.0, 1.02),fancybox=True, fontsize=10, labelspacing=0.1)

plt.title('MJO LPT System Tracks: 1998 - 2018\n' +"(N = " + str(grand_total_count) + " MJO, green; LPTs, gray)" )



##
## Plot track density
##
ax2 = fig.add_subplot(212)


# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map2 = Basemap(llcrnrlon=0.,llcrnrlat=-50.,urcrnrlon=360.,urcrnrlat=50.,\
              resolution='l',projection='merc')
# draw coastlines, country boundaries, fill continents.
map2.fillcontinents(color='wheat',zorder=1)
map2.drawcoastlines(color='k',linewidth=0.7,zorder=10)
map2.drawcountries(color='k',linewidth=0.5,zorder=10)
# draw the edge of the map projection region (the projection limb)
map2.drawmapboundary()
# draw lat/lon grid lines every 30 degrees.
map2.drawmeridians(np.arange(0,360,30),zorder=15)
map2.drawparallels(np.arange(-90,90,30),zorder=15)


ncfile='../../data/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/postproc/lpt_track_density.nc'

DS = Dataset(ncfile,'r')

lon = DS['lon'][:]
lat = DS['lat'][:]
density = DS['track_density_mjo_eastward_portion'][:]
x1, y1 = np.meshgrid(lon, lat)
x2, y2 = map2(x1+1.125, y1+1.125)

plt.contour(x2, y2, density, np.arange(5,75,5), cmap="jet")

plt.title('MJO Track Density: 1998 - 2018')

xt, yt = map2(160, 41)
plt.text(xt, yt, 'Contours every 5\n 2.5 deg. boxes', fontsize=10)


plt.tight_layout()

    
print(fout)
plt.savefig(fout, dpi = 150,bbox_inches='tight')


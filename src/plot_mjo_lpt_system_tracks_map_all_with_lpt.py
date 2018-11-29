import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import glob
import os

dir = '../data/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/timeclusters'
out_dir = '../plots/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/systems'
plt.close('all')


## Read in MJO LPTs
"""
mc_crossing_file = '/home/orca/bkerns/data/lpt/trmm/processed/g20_72h/thresh12/identify_eastward_propagation/list_mc_crossing_io_lpts.txt'
non_mc_crossing_file = '/home/orca/bkerns/data/lpt/trmm/processed/g20_72h/thresh12/identify_eastward_propagation/list_non_mc_crossing_io_lpts.txt'
wpac_file = '/home/orca/bkerns/data/lpt/trmm/processed/g20_72h/thresh12/identify_eastward_propagation/list_wpac_lpts.txt'
"""
mjo_file = '/home/orca/bkerns/data/lpt/trmm2/processed/g20_72h/thresh12/identify_eastward_propagation/mjo_lpt_list.rejoin2.txt'
clumps_file = '/home/orca/bkerns/data/lpt/trmm2/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'

D_mjo = np.loadtxt(mjo_file, skiprows=1)

#D_mc_crossing = np.loadtxt(mc_crossing_file, skiprows=1)
#D_non_mc_crossing = np.loadtxt(non_mc_crossing_file, skiprows=1)
#D_wpac = np.loadtxt(wpac_file, skiprows=1)
D_all = D_mjo   #np.vstack((D_mc_crossing, D_non_mc_crossing, D_wpac))
mjo={}
mjo['year'] = D_all[:,0]
mjo['idx'] = D_all[:,1]

D_clumps = np.loadtxt(clumps_file, skiprows=1)
clumps={}
clumps['year'] = D_clumps[:,0]
clumps['idx'] = D_clumps[:,1]
clumps['clump_num'] = D_clumps[:,2]

os.system('mkdir -p ' + out_dir)


grand_total_count = 0
grand_total_lpt_count = 0



## Set up figure
    
fig = plt.figure(figsize=(11,8))
ax1 = fig.add_subplot(111)
    
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


fout = (out_dir + '/lpt_mjo_system_tracks_20yrs__map_with_lpt.png')


for fn in sorted(glob.glob(dir+'/TIMECLUSTERS_*.rejoin2.nc')):

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

        #print(cnum)

        ## mask out non-LPT systems.
        for idx in np.unique(clumps_this_year['idx'][clumps_this_year['clump_num'] == cnum]):

            #print(idx)

            lon_this = lon[lptid == idx]
            lat_this = lat[lptid == idx]
            
            if np.logical_and(mjo['year'] == year, mjo['idx'] == idx).any():


                mjo_count += 1
                
                x, y = map(lon_this, lat_this)
                
                map.plot(x,y,'k',zorder=5)
                map.plot(x[0], y[0], 'ro', markersize=3, zorder=1000)
                map.plot(x[-1], y[-1], 'mx', markersize=5, zorder=1000, markeredgewidth=2.0)
                
                break  # Just take the first one for now. TODO: This needs to be the one with longest eastward propagation period.

            ## If I get to this point, not an MJO LPT.
            x, y = map(lon_this, lat_this)
            
            map.plot(x,y,zorder=2, color=[0.5, 0.5, 0.5], linewidth=0.5)
            
            
    print(mjo_count)
    grand_total_count += mjo_count

    
xt, yt = map(170.0, 40.0)

plt.text(xt, yt, ("N = " + str(grand_total_count) + "MJO, black\n" + "(of " + str(grand_total_lpt_count) + "LPTs, gray)")  , fontsize=10, fontweight='bold')


plt.title('MJO LPT System Tracks: 1998 - 2018')

print(fout)
plt.savefig(fout, dpi = 150,bbox_inches='tight')




    

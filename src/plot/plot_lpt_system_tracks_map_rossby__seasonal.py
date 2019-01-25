import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import glob
import os
import datetime as dt

#months_to_keep = [12,1,2]
#season_label = "djf"
months_to_keep = [6,7,8]
season_label = "jja"

##########################################################

dir = '../../data/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/timeclusters'
out_dir = '../../plots/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/systems'
plt.close('all')


## Read in MJO LPTs
mjo_file = '/home/orca/bkerns/data/lpt/trmm3/processed/g20_72h/thresh12/identify_eastward_propagation/rossby_lpt_list.rejoin2.txt'
clumps_file = '/home/orca/bkerns/data/lpt/trmm3/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'

D_rossby = np.loadtxt(mjo_file, skiprows=1)
rossby={}
rossby['year'] = D_rossby[:,0]
rossby['idx'] = D_rossby[:,1]
rossby['start'] = D_rossby[:,8]
rossby['end'] = D_rossby[:,9]


D_clumps = np.loadtxt(clumps_file, skiprows=1)
clumps={}
clumps['year'] = D_clumps[:,0]
clumps['idx'] = D_clumps[:,1]
clumps['clump_num'] = D_clumps[:,2]

os.system('mkdir -p ' + out_dir)


grand_total_rossby_count = 0
grand_total_lpt_count = 0



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


fout = (out_dir + '/rossby_system_tracks_20yrs__map__'+season_label+'.png')


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
    stitchtime = DS['end_of_accumulation_time'][:]

    rossby_count = 0
    lpt_count = 0
    
    for cnum in sorted(np.unique(clumps_this_year['clump_num'])):

        lpt_count += 1
        grand_total_lpt_count += 1

        ## mask out non-LPT systems.
        for idx in np.unique(clumps_this_year['idx'][clumps_this_year['clump_num'] == cnum]):


            lon_this = lon[lptid == idx]
            lat_this = lat[lptid == idx]
            
            ## Filter by months
            ## Skip this LPT if it does not have any times within the "months_to_keep"
            stitchtime_this = stitchtime[lptid == idx]
            dt_this = [dt.datetime(1970,1,1,0,0,0) + dt.timedelta(hours=x) for x in stitchtime_this]
            months_this = [x.month for x in dt_this]
            if len(set(months_this).intersection(months_to_keep)) < 1:
                continue


            ## grey track in background.
            x, y = map(lon_this, lat_this)            
            #h00 = map.plot(x,y,zorder=2, color=[0.5, 0.5, 0.5], linewidth=0.7)
            
            if np.logical_and(rossby['year'] == year, rossby['idx'] == idx).any():

                rossby_count += 1
                grand_total_rossby_count += 1

                idx11 = int(rossby['start'][np.logical_and(rossby['year'] == year, rossby['idx'] == idx)][0])
                idx22 = int(rossby['end'][np.logical_and(rossby['year'] == year, rossby['idx'] == idx)][0])
                
                x, y = map(lon_this[idx11:idx22], lat_this[idx11:idx22])
                                
                h0=map.plot(x, y, color='yellowgreen', zorder=5, linewidth=0.7)
                h1=map.plot(x[0], y[0], 'ko', markersize=2, zorder=1000)
                h2=map.plot(x[-1], y[-1], 'rx', markersize=5, zorder=1000, markeredgewidth=1.5)                            
            
    print(rossby_count)


leg = plt.legend((h0[0], h1[0], h2[0],),('Westward LPT', 'Westward Start','Westward End',), loc=(0.0, 1.02),fancybox=True, fontsize=10, labelspacing=0.1)

plt.title('Westward LPT System Tracks: 1998 - 2018 ('+season_label.upper()+')\n' + "(N = " + str(grand_total_rossby_count) + ")")



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


ncfile='../../data/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/postproc/lpt_track_density__'+season_label+'.nc'

DS = Dataset(ncfile,'r')

lon = DS['lon'][:]
lat = DS['lat'][:]
density = DS['track_density_rossby_westward_portion'][:]
x1, y1 = np.meshgrid(lon, lat)
x2, y2 = map2(x1+1.125, y1+1.125)

plt.contour(x2, y2, density, np.arange(2,72,2), cmap="jet")

plt.title('Westward LPTs Track Density: 1998 - 2018 ('+season_label.upper()+')')

xt, yt = map2(160, 41)
plt.text(xt, yt, 'Contours every 2\n 2.5 deg. boxes', fontsize=10)

plt.tight_layout()






print(fout)
plt.savefig(fout, dpi = 150,bbox_inches='tight')




    

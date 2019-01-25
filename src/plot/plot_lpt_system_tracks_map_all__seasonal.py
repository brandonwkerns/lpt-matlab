import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import glob
import os
import datetime as dt

months_to_keep = [12,1,2]
season_label = "djf"
#months_to_keep = [6,7,8]
#season_label = "jja"

##########################################################

dir = '../../data/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/timeclusters'
out_dir = '../../plots/trmm_keep_overlapping_tracks_no_jump_cutoff/processed/g20_72h/thresh12/systems'
plt.close('all')
clumps_file = '/home/orca/bkerns/data/lpt/trmm3/processed/g20_72h/thresh12/identify_eastward_propagation/clumps_of_worms.rejoin2.txt'

D_clumps = np.loadtxt(clumps_file, skiprows=1)
clumps={}
clumps['year'] = D_clumps[:,0]
clumps['idx'] = D_clumps[:,1]
clumps['clump_num'] = D_clumps[:,2]

os.system('mkdir -p ' + out_dir)


grand_total_count = 0
grand_total_lpt_count = 0
grand_total_clumps_count = 0


fout = (out_dir + '/lpt_system_tracks_20yrs__map__'+season_label+'.png')

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
    stitchtime = DS['end_of_accumulation_time'][:]

    clumps_count = 0
    lpt_count = 0
    
    for cnum in sorted(np.unique(clumps_this_year['clump_num'])):

        grand_total_clumps_count += 1

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

            lpt_count += 1
            grand_total_lpt_count += 1
            
            x, y = map(lon_this, lat_this)
            
            h0=map.plot(x,y,zorder=2, color=[0.5, 0.5, 0.5], linewidth=0.7)
            h1=map.plot(x[0],y[0],'ko',zorder=10,  markersize=2, linewidth=0.5)
            h2=map.plot(x[-1],y[-1],'rx',zorder=10,  markersize=5.0, markeredgewidth=1.5)
 

leg = plt.legend((h0[0], h1[0], h2[0],),('LPT', 'Initiation','Dissipation',), loc=(0.0, 1.02),fancybox=True, fontsize=10, labelspacing=0.1)


plt.title('LPT System Tracks: 1998 - 2018 ('+season_label.upper()+')\n' + "(" + str(grand_total_lpt_count) + " LPTs)")


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
density = DS['track_density_lpt'][:]
x1, y1 = np.meshgrid(lon, lat)
x2, y2 = map2(x1, y1)

plt.contour(x2, y2, density, np.arange(5,75,5), cmap="jet")

plt.title('LPT Track Density: 1998 - 2018 ('+season_label.upper()+')')

xt, yt = map2(160, 41)
plt.text(xt, yt, 'Contours every 5\n 2.5 deg. boxes', fontsize=10)

plt.tight_layout()


print(fout)
plt.savefig(fout, dpi = 150,bbox_inches='tight')




    

import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import glob
import os

dir = '../data/trmm/processed/g20_72h/thresh12/objects'
out_dir = '../plots/trmm/processed/g20_72h/thresh12/objects'
plt.close('all')

os.system('mkdir -p ' + out_dir)

for fn in sorted(glob.glob(dir+'/objects_*.nc')):

    print(fn, flush=True)

    out_str = fn[-24:-3]
    fout = (out_dir + '/lp_objects_' + out_str + '__map.png')

    ## Set up figure
    fig = plt.figure(figsize=(11,8))
    ax1 = fig.add_subplot(111)



    ## Read
    DS = Dataset(fn)

    lon = DS['lon'][:]
    lat = DS['lat'][:]
    time = DS['time'][:]
    area = DS['area'][:] / (1e5)

    ## Plot
    # set up orthographic map projection with
    # perspective of satellite looking down at 50N, 100W.
    # use low resolution coastlines.
    #map = Basemap(projection='hammer',lon_0=180)
    map = Basemap(llcrnrlon=0.,llcrnrlat=-50.,urcrnrlon=360.,urcrnrlat=50.,\
                resolution='l',projection='merc')
    # draw coastlines, country boundaries, fill continents.
    map.fillcontinents(color='lightgray',zorder=1)
    map.drawcoastlines(color='k',linewidth=0.7,zorder=10)
    map.drawcountries(color='k',linewidth=0.5,zorder=10)
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary()
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30),zorder=15)
    map.drawparallels(np.arange(-90,90,30),zorder=15)


    x, y = map(lon, lat)
    map.scatter(x,y, area, marker='o',facecolors='none', edgecolors='blue',zorder=5)

    plt.title(out_str)

    print(fout)
    plt.savefig(fout, dpi = 150,bbox_inches='tight')

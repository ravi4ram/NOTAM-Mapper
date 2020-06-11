#  notam-mapper.py
#  
#  Copyright 2020 ravi <ravi@CE24X005>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

#  returns [lat,lon] of launch pad coord and danger zone coords
def notam_coordinate_parser(notam_file_name):
	num_lines = sum(1 for line in open(notam_file_name))
	
	with open(notam_file_name,"r") as fi:
		# fsm flags
		flag_launchpad = False; flag_dng = False; 
		flag_coordinates = False; flag_status = False;
		# launch pad coordinates
		lat_pad = 0; lon_pad = 0;
		# counters for danger zones
		counter_dng = 0; counter_dng_prev = 0;
		# list of lats and longs of each danger zone
		lat = 0; lon = 0; lats=[]; longs=[];
		# list of [lats,longs]
		polygons=[]; 
		
		for i, ln in enumerate(fi, start=1):
			# strip leading spaces & trailing junks
			ln = ln.strip(); ln = ln.strip('\n'); ln = ln.strip('\t')

			key_start = ["LAUNCH PAD COORD"]
			if any(c in ln for c in key_start):							
				flag_launchpad = True
			
			key_danger = ["DNG ZONE", "DANGER ZONE"]
			if any(c in ln for c in key_danger):
				flag_dng = True; flag_coordinates = False; 
				# increment dng zone counter
				counter_dng = counter_dng + 1
			
			# set termination condition
			key_termination = ["ALL COORD ARE IN DEG AND MIN", "CLOSURE/ALTERNATE"]
			if any(c in ln for c in key_termination):			
				# include the final parsed list before exiting
				# append to polygen list(if not empty) 
				if (lats and longs):
					polygons.append([lats,longs])
				# reset the flags
				flag_status = True; flag_dng = False; 
				flag_coordinates = False; flag_launchpad = False;
			# finished extraction and break out
			if(flag_status): break
			# extract coordinates
			if(flag_launchpad):
				flag_launchpad = False;
				# with decimal part
				pat  = r'[0-9]{4,7}\.?\d+[NSEW]' 
				lst = re.findall(pat, ln)
				if(len(lst) > 0):
					lat_pad, lon_pad = notam_dms2dd(lst[0]), notam_dms2dd(lst[1])	
								
			if(flag_dng):
				pat = r'[0-9]{4,7}[NSEW]'
				lst = re.findall(pat, ln)
				# in case if dzone coordinates are in multiple lines 
				if(len(lst) == 2):
					lat, lon = notam_dms2dd(lst[0]), notam_dms2dd(lst[1])

					if( counter_dng_prev == counter_dng ):
						lats.append(lat); longs.append(lon);
					else:
						counter_dng_prev = counter_dng
						# append to polygen list(if not empty) 
						if (lats and longs):
							polygons.append([lats,longs])
						# and clear the old ones
						lats = []; longs = [];
						lats.append(lat); longs.append(lon);	
				# in case if dzone coordinates are in the same line 
				if(len(lst) == 8):
					for i in range(0,8,2):
						lat, lon = notam_dms2dd(lst[i]), notam_dms2dd(lst[i+1])

						if( counter_dng_prev == counter_dng ):
							lats.append(lat); longs.append(lon);
						else:
							counter_dng_prev = counter_dng
							# append to polygen list(if not empty) 
							if (lats and longs):
								polygons.append([lats,longs])
							# and clear the old ones
							lats = []; longs = [];
							lats.append(lat); longs.append(lon);	
			# append to polygen list(file last line-no termination condition) 
			if (num_lines==i and lats and longs):
				polygons.append([lats,longs])				
	# done
	return [lat_pad, lon_pad], polygons;

#  returns lat/lon in decimal (ddmmss to decimal conversion)
def notam_dms2dd(s):
	# direction
	direction = s[-1]	
	# remove direction
	s = s[:-1]
	# get string length
	size = len(s)	

	if (size % 2) == 0:	
		degrees = s[0:2]; minutes = s[2:4]; seconds = s[4:6];			
	else:
		degrees = s[0:3]; minutes = s[3:5]; seconds = s[5:7];
	
	if( size == 10):
		degrees = s[0:3]; minutes = s[3:5]; seconds = s[5:10];		
	if( size == 9):
		degrees = s[0:2]; minutes = s[2:4]; seconds = s[4:9];	

	if( (size == 4) or (size == 5) ): seconds = 0.0

	# compute decimal coordinate with max 5 decimals
	dd = float("{0:.5f}".format( float(degrees) + float(minutes)/60 + float(seconds)/(60*60) ) );
	if direction == 'S' or direction == 'W': dd *= -1

	return dd;

#  draw a merc map with given boundaries
def drawBasemap(northbounds, eastbounds, southbounds, westbounds):
	fig = plt.figure()
	mp = Basemap( projection='merc', resolution=None, \
	area_thresh = 0.1, ax=fig.gca(), \
	urcrnrlat=northbounds, urcrnrlon=eastbounds, \
	llcrnrlat=southbounds, llcrnrlon=westbounds )
	mp.bluemarble(scale=0.98)
	mp.drawparallels(np.arange(-90.,91.,7.), color='black', labels=[True,True,False,False],dashes=[2,2])
	mp.drawmeridians(np.arange(-180.,181.,10.),color='black', labels=[False,False,False,True],dashes=[2,2])
	return mp

#  draws colored polygon-danger zones-at the given lat,lon
def map_danger_zones( m, lats, lons, mpoly ):
	launch_path_lat=[]; launch_path_lon=[];
	# launch pad
	launch_path_lat.append(lat)
	launch_path_lon.append(lon)
	for i, (lats, lons) in enumerate(mpoly, start=2):
		x, y = m( lons, lats )
		xy = zip(x,y)
		poly = Polygon( list(xy), facecolor='red', alpha=0.4 )
		plt.gca().add_patch(poly)	
		plt.text(x[0]+20, y[0], '  D'+str(i), color='white', fontsize=7)		
		# center points
		lat_c, lon_c = get_center_dzones(lats, lons)
		launch_path_lat.append(lat_c)
		launch_path_lon.append(lon_c)		
	# connect the paths
	x, y = map(launch_path_lon, launch_path_lat)
	m.plot(x, y, linestyle=':', color='m')
	return

#  draws a circle danger zone around launchpad
def map_launch_pad(lat_0, lon_0, m, col='r'):	
	x, y = m(lon_0, lat_0)
	plt.text(x, y, '   Sriharikota', color='white', fontsize=7)

	dist = 10
	dist = dist * 1852  
	theta_arr = np.linspace(0, np.deg2rad(360), 100)
	lat_0 = np.deg2rad(lat_0)
	lon_0 = np.deg2rad(lon_0)

	coords_new = []
	for theta in theta_arr:
		coords_new.append(calc_new_coord(lat_0, lon_0, theta, dist))
	lat = [item[0] for item in coords_new]
	lon = [item[1] for item in coords_new]

	x, y = m(lon, lat)
	m.plot(x, y, col)
	return

"""
Calculate coordinate pair given starting point, radial and distance
Method from: http://www.geomidpoint.com/destination/calculation.html
"""
def calc_new_coord(lat1, lon1, rad, dist):
    flat = 298.257223563
    a = 2 * 6378137.00
    b = 2 * 6356752.3142

    # Calculate the destination point using Vincenty's formula
    f = 1 / flat
    sb = np.sin(rad)
    cb = np.cos(rad)
    tu1 = (1 - f) * np.tan(lat1)
    cu1 = 1 / np.sqrt((1 + tu1*tu1))
    su1 = tu1 * cu1
    s2 = np.arctan2(tu1, cb)
    sa = cu1 * sb
    csa = 1 - sa * sa
    us = csa * (a * a - b * b) / (b * b)
    A = 1 + us / 16384 * (4096 + us * (-768 + us * (320 - 175 * us)))
    B = us / 1024 * (256 + us * (-128 + us * (74 - 47 * us)))
    s1 = dist / (b * A)
    s1p = 2 * np.pi

    while (abs(s1 - s1p) > 1e-12):
        cs1m = np.cos(2 * s2 + s1)
        ss1 = np.sin(s1)
        cs1 = np.cos(s1)
        ds1 = B * ss1 * (cs1m + B / 4 * (cs1 * (- 1 + 2 * cs1m * cs1m) - B / 6 * \
            cs1m * (- 3 + 4 * ss1 * ss1) * (-3 + 4 * cs1m * cs1m)))
        s1p = s1
        s1 = dist / (b * A) + ds1

    t = su1 * ss1 - cu1 * cs1 * cb
    lat2 = np.arctan2(su1 * cs1 + cu1 * ss1 * cb, (1 - f) * np.sqrt(sa * sa + t * t))
    l2 = np.arctan2(ss1 * sb, cu1 * cs1 - su1 * ss1 * cb)
    c = f / 16 * csa * (4 + f * (4 - 3 * csa))
    l = l2 - (1 - c) * f * sa * (s1 + c * ss1 * (cs1m + c * cs1 * (-1 + 2 * cs1m * cs1m)))
    d = np.arctan2(sa, -t)
    finaltc = d + 2 * np.pi
    backtc = d + np.pi
    lon2 = lon1 + l
    return (np.rad2deg(lat2), np.rad2deg(lon2))

def get_center_dzones(lat1,lon1):    
	# input in degrees
	if (len(lat1) <= 0):
		return false;

	num_coords = len(lat1)

	X = 0.0; Y = 0.0; Z = 0.0

	for i in range (len(lat1)):
		lat = lat1[i] * np.pi / 180
		lon = lon1[i] * np.pi / 180

		a = np.cos(lat) * np.cos(lon)
		b = np.cos(lat) * np.sin(lon)
		c = np.sin(lat);

		X += a
		Y += b
		Z += c

	X /= num_coords
	Y /= num_coords
	Z /= num_coords

	lon = np.arctan2(Y, X)
	hyp = np.sqrt(X * X + Y * Y)
	lat = np.arctan2(Z, hyp)

	newX = (lat * 180 / np.pi)
	newY = (lon * 180 / np.pi)
	return newX, newY	


# __main method__ 
if __name__=="__main__": 
	# input notam file name 
	notam_filename = "notam-gslv-m1.txt" 
	 
	pos, poly = notam_coordinate_parser(notam_filename)
	lat = pos[0]; lon = pos[1]; 
	# map boundaries
	llcrnrlat=-15; urcrnrlat=20; llcrnrlon=70; urcrnrlon=110;
	map = drawBasemap(urcrnrlat, urcrnrlon, llcrnrlat, llcrnrlon)
	map_launch_pad(lat, lon, map, 'r')
	map_danger_zones(map, lat, lon, poly)
	# show it
	plt.title(Path(notam_filename).stem)
	plt.savefig(Path(notam_filename).stem+'.png')
	plt.show()

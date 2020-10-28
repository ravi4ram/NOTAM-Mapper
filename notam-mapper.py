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
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Polygon
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def notam_coordinate_parser(notam_file_name):
	# number of lines
	num_lines = sum(1 for line in open(notam_file_name))
	
	with open(notam_file_name,"r") as fi:
		# fsm flags
		flag_launchpad = False; flag_dng = False; 
		flag_coordinates = False; flag_status = False;
		# launch pad coordinates
		lat_pad = 0; lon_pad = 0;
		# counters for danger zones
		counter_dng = 0;
		counter_dng_prev = 0;
		# list of lats and longs of each danger zone
		lat = 0; lon = 0;
		lats=[]; longs=[];
		# list of [lats,longs]
		polygons=[]; 
		
		for i, ln in enumerate(fi, start=1):
			# strip leading spaces
			ln = ln.strip()
			# strip trailing junks
			ln = ln.strip('\n'); ln = ln.strip('\t')
							
			if ln.startswith("LAUNCH PAD COORD"):
				flag_launchpad = True
			
			key_danger = ["DNG ZONE", "DANGER ZONE"]
			if any(c in ln for c in key_danger):
			# if ln.startswith( ("DNG ZONE", "DANGER ZONE") ):	
				flag_dng = True; flag_coordinates = False; 
				# increment dng zone counter
				counter_dng = counter_dng + 1
			
			# set termination condition
			key_termination = ["ALL COORD ARE IN DEG AND MIN", "CLOSURE/ALTERNATE"]
			if any(c in ln for c in key_termination):			
			# if ln.startswith( ("ALL COORD ARE IN DEG AND MIN", "CLOSURE/ALTERNATE") ):
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
				# print(ln)
				flag_launchpad = False;
				# pat  = r'[0-9]{4,7}[NSEW]'
				# with decimal part
				pat  = r'[0-9]{4,7}\.?\d+[NSEW]' 
				lst = re.findall(pat, ln)
				if(len(lst) > 0):
					lat_pad, lon_pad = notam_dms2dd(lst[0]), notam_dms2dd(lst[1])	
								
			if(flag_dng):
				pat  = r'[0-9]{4,7}[NSEW]'
				# pat  = r'[0-9]{4,7}\.?\d+[NSEW]' 
				lst = re.findall(pat, ln)
				# case if dzone coordinates are in multiple lines 
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
				# case if dzone coordinates are in the same line 
				if(len(lst) == 8):
					for i in range(0,8,2):
						lat, lon = notam_dms2dd(lst[i]), notam_dms2dd(lst[i+1])
						# print(lat,"------",lon)

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
	
def notam_dms2dd(s):
	# copy last digit as direction
	direction = s[-1]	
	# remove direction
	s = s[:-1]
	# get string length
	size = len(s)	
	# check if len is odd/even.  odd will have 3-digit degrees
	if (size % 2) == 0:	
		degrees = s[0:2]; minutes = s[2:4]; seconds = s[4:6];			
	else:
		degrees = s[0:3]; minutes = s[3:5]; seconds = s[5:7];

	# special cases with seconds in float	
	if( size == 10):
		degrees = s[0:3]; minutes = s[3:5]; seconds = s[5:10];		
	if( size == 9):
		degrees = s[0:2]; minutes = s[2:4]; seconds = s[4:9];	

	if( (size == 4) or (size == 5) ): seconds = 0.0

	# compute decimal coordinate with max 5 decimals
	dd = float("{0:.5f}".format( float(degrees) + float(minutes)/60 + float(seconds)/(60*60) ) );
	if direction == 'S' or direction == 'W': dd *= -1

	return dd;

# polygen centroid
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

#  draw a PlateCarree map with given boundaries
def map_launch_area(lllon, urlon, lllat, urlat):
	# center the map on
	c_longitude = 80
	# subplt
	fig, ax = plt.subplots(figsize=(8, 6)) #
	fig.subplots_adjust(bottom=0.01)
	fig.tight_layout()	
	# projection	
	ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=c_longitude))
	ax.set_global()
	# zoom in
	ax.set_extent([lllon, urlon, lllat, urlat], crs=ccrs.PlateCarree())
	# gridlines
	gl1 = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=1, linestyle='--', edgecolor='dimgrey')	
	gl1.top_labels = False; gl1.left_labels = True
	# add coastlines for reference                                                                                                
	ax.coastlines(resolution='50m')		
	ax.add_feature(cfeature.OCEAN)
	ax.stock_img()
	return ax

#  draws colored polygon-danger zones-at the given lat,lon
def map_danger_zones( ax, lat, lon, mpoly ):
	launch_path_lat=[]; launch_path_lon=[];
	# launch pad
	ax.scatter( lon, lat, transform=ccrs.PlateCarree(), s=20, color='red')
	ax.text(lon, lat, '   Sriharikota', color='k', fontsize='x-small', transform=ccrs.PlateCarree()) 	
	# path co-ordinates
	launch_path_lat.append(lat)
	launch_path_lon.append(lon)
	# danger zone areas
	for i, (lats, lons) in enumerate(mpoly, start=2):
		xy = zip(lons, lats )
		poly = Polygon( list(xy), transform=ccrs.PlateCarree(), closed=True, fill=True, fc='red', ec='m', lw=2, alpha=0.4 )
		ax.add_patch(poly)	
		ax.text( lons[0], lats[0], '  D'+str(i), color='k', fontsize='x-small', transform=ccrs.PlateCarree()) 	
		# center points
		lat_c, lon_c = get_center_dzones(lats, lons)
		launch_path_lat.append(lat_c)
		launch_path_lon.append(lon_c)		
	# connect the paths
	ax.plot(launch_path_lon, launch_path_lat, transform=ccrs.PlateCarree(), linestyle=':', color='b')
	return	

# __main method__ 
if __name__=="__main__": 
	# input notam file name 
	notam_filename = "./test/notam-gslv-f10.txt"  
	# parser returns launch pad coord and danger zone coords.	 
	pos, poly = notam_coordinate_parser(notam_filename)
	# launchpad co-ordinates
	lat = pos[0]; lon = pos[1];
	# map boundaries
	llcrnrlon=70; urcrnrlon=110; llcrnrlat=-15; urcrnrlat=20; 
	# draw map
	ax = map_launch_area(llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat)
	# draw danger zones
	map_danger_zones(ax, lat, lon, poly)
	# show it
	plt.title(Path(notam_filename).stem)
	plt.savefig(Path(notam_filename).stem+'.png')
	plt.show()

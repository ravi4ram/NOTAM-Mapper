# NOTAM-Mapper
> This script reads the NOTAM (Notice To Airmen) file and generates a map with danger zones marked on it. 

## Table of contents
* [General info](#general-info)
* [Screenshots](#screenshots)
* [Setup](#setup)
* [How to run ](#how)
* [To-do list](#to-do)

## General info
Current version of this project is written for rocket launches from ISRO (Indian Space Research Organization). NOTAM'S issued for few ISRO launches were used as test cases and are included under the directory named 'test' 

## Screenshots
Outputs produced by the script in running with the NOTAMs issued for the launches of PSLV-C47, PSLV-C48, GSLV-F10 and GSLV-M1.    

![results](./img/results_merged.png)   


## Setup
Script is written with python (Version: 3.6) on linux. Additional modules required:   

* numpy  (tested with Version: 1.18.4 )
* matplotlib  (tested with Version: 2.1.1 )
* basemap ( tested with Version: 1.1.0 )

## How to run   
* copy one of the notams into a txt file in the script directory, say notam-LV-noXX.txt
* edit code at line no.286 (around) to read the new NOTAM.  
 `notam_filename = "notam-LV-noXX.txt"`   
* run `python notam-mapper.py`. It
generates a png file and shows the plot.

## To-do list
* This initial version is ISRO specific and the map boundaries and hard-coded. This could be made automatic by taking min-max values of latitude/longitude read through the danger zones. Another way is to store specific boundaries for each launch port in a dictionary and display them accordingly. 
* Need to test with NOTAMS of others launch agencies and modify the code accordingly.    
* Script could be modified to be as a bot running for twitter, reddit etc.


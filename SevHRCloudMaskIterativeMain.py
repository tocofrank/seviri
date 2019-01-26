# ========================================================
# Import modules
# ========================================================

#Import python numerical
import numpy as np

#Import python data structures and data statistics tools
import pandas as pd

#Import python scientific tools
from scipy.signal import tukey
from scipy.interpolate import RegularGridInterpolator as gridintp

#Import python H5
import h5py

#Import system and os tools
import sys
import os
from pathlib2 import Path

#Import SEVIRI subroutines
sys.path.append('../subroutines/')
from SeviriSubroutines import *

#Import basic functions
from Functions import *








# ========================================================
# Define directories, the year, x-day average for the clear sky composite, start and end time, time interval (in minutes) for the temporal resolution
# ========================================================

l1b_inpath = "/vols/altair/datasets/eumcst/msevi_rss/l15_hdf/eu/"
ct_inpath = "/vols/talos/home/fabian/proj/2016-06_hdcp2_omodule/nwcsaf_products/"
csc_inpath = "/vols/talos/home/frank/data/seviri/clear_sky_composites/"
thresh_outpath = "/vols/talos/home/frank/data/seviri/hrv_cm/"

date_year = 2013
date_month_start = 4 # usually 1
date_day_start = 1
julian_start = JulianDay(date_year, date_month_start, date_day_start)
date_month_end = 7 # usually 12
date_day_end = 31
julian_end = JulianDay(date_year, date_month_end, date_day_end)

julian_step = 16

th_min = 8
tm_min = 0
th_max = 16
tm_max = 0

dt = 30

cscth_min = 8
csctm_min = 0
cscth_max = 16
csctm_max = 0

cscdt = 60

# Generate the array, which specifies which Julian days are included in each clear sky composite file; the last date is set manually, depending on julian_end, thus this last bin is usually a little larger
julian_array = np.arange(julian_start, (np.round(julian_end/julian_step)+1)*julian_step, julian_step, dtype=int)

if julian_array[len(julian_array)-1] < julian_end:
    julian_array[len(julian_array)-1] = julian_end

# Generate the time array, which specifies the time steps for each clear sky composite
time_array = np.zeros( (np.int((th_max + tm_max/60. - (th_min+tm_min/60.))/(dt/60.))+1,2), dtype=np.int )

for i in range(0, len(time_array)):
    time_array[i,0] = np.int(th_min + i*dt/60.)
    time_array[i,1] = np.int( np.round((th_min + i*dt/60. - time_array[i,0])*60.) )

csctime_array = np.zeros( (np.int((cscth_max + csctm_max/60. - (cscth_min+csctm_min/60.))/(cscdt/60.))+1,2), dtype=np.int )

for i in range(0, len(csctime_array)):
    csctime_array[i,0] = np.int(cscth_min + i*cscdt/60.)
    csctime_array[i,1] = np.int( np.round((cscth_min + i*cscdt/60. - csctime_array[i,0])*60.) )








# ========================================================
# Loop over all elements of julian_array (this is the number of files we will get) and time_array (these are the time steps in each file)
# ========================================================

for ff in range(0, len(julian_array)-1):

	area = "EU" # Europe
	dimy = 1800
	dimx = 2400
	area = "DE" # Germany
	dimy = 720
	dimx = 1200

	sev_dailyhrv = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.float32)
	sev_dailyct = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.int)
	sev_dailyctc = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.int)
	sev_dailyoricm = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.int)

	landthreshold = np.zeros((1,len(csctime_array)), dtype=np.float32)
	seathreshold = np.zeros((1,len(csctime_array)), dtype=np.float32)
	coastthreshold = np.zeros((1,len(csctime_array)), dtype=np.float32)

	# Read clear sky composite file
	if ff < len(julian_array)-2:
		sevclearfilename_dummy11, sevclearfilename_dummy12 = JulianDayToDate(date_year, julian_array[ff])
		sevclearfilename_dummy21, sevclearfilename_dummy22 = JulianDayToDate(date_year, julian_array[ff+1])
	else:
		sevclearfilename_dummy11, sevclearfilename_dummy12 = JulianDayToDate(date_year, julian_array[ff])
		sevclearfilename_dummy21, sevclearfilename_dummy22 = JulianDayToDate(date_year, julian_array[len(julian_array)-1])
	if sevclearfilename_dummy11 < 10:
		sevclearfilename_dummy11 = "0"+str(sevclearfilename_dummy11)
	else:
		sevclearfilename_dummy11 = str(sevclearfilename_dummy11)
	if sevclearfilename_dummy12 < 10:
		sevclearfilename_dummy12 = "0"+str(sevclearfilename_dummy12)
	else:
		sevclearfilename_dummy12 = str(sevclearfilename_dummy12)
	if sevclearfilename_dummy21 < 10:
		sevclearfilename_dummy21 = "0"+str(sevclearfilename_dummy21)
	else:
		sevclearfilename_dummy21 = str(sevclearfilename_dummy21)
	if sevclearfilename_dummy22 < 10:
		sevclearfilename_dummy22 = "0"+str(sevclearfilename_dummy22)
	else:
		sevclearfilename_dummy22 = str(sevclearfilename_dummy22)

	sevclearfilename = csc_inpath + "msg2-sevi-"+str(date_year)+"."+sevclearfilename_dummy11+sevclearfilename_dummy12+"-"+sevclearfilename_dummy21+sevclearfilename_dummy22+"-clearskycomp-rss-eu.c2.h5"
	csc_data = SevReadClearSkyCompfile(sevclearfilename, csctime_array)
	
	# Reduce dataset to Germany domain if wanted
	if area == "DE":
		y0 = 60
		y1 = 300
		x0 = 320
		x1 = 720
				
		sy = slice(y0*3,y1*3)
		sx = slice(x0*3,x1*3)
		csc_data = np.array(csc_data[sy,sx,:])
		
	# ========================================================
	# Iterative loop to adjust the clear sky composite
	# ========================================================

	it_loop_land = 10
	it_loop_sea = 10
	it_loop_coast = 10
	
	landthreshold_old = 10
	seathreshold_old = 10
	coastthreshold_old = 10

	for it_loop in range(0, 10):
	
		if it_loop > 0:
			landthreshold_old = landthreshold
			seathreshold_old = seathreshold
			coastthreshold_old = coastthreshold
	
		if it_loop_land > 0.01 and it_loop_sea > 0.01 and it_loop_coast > 0.01:

			# Loop over all time steps
			ttt = 0 # this is the index for the csctime_array (i.e., 8 am, 9am, ...)
			tt_index = -1 # this is the counter for the time_array (i.e., 8:00am, 8:05 am, 8:10 am), which is 60./np.float(dt) long for each element of csctime_array (except for the first and last, which are half that size (e.g., 8:00-8:25)

			tt = 0
			while tt < len(time_array):
	
				tt_index += 1

				for fff in range(0, julian_step):
	
					date_month, date_day = JulianDayToDate(date_year, julian_array[ff]+fff)
					date_hour = time_array[tt,0]
					date_minute = time_array[tt,1]
					date_year = str(date_year)
					if date_month < 10:
						date_month = "0"+str(date_month)
					else:
						date_month = str(date_month)
					if date_day < 10:
						date_day = "0"+str(date_day)
					else:
						date_day = str(date_day)
					if date_hour < 10:
						date_hour = "0"+str(date_hour)
					else:
						date_hour = str(date_hour)
					if date_minute < 10:
						date_minute = "0"+str(date_minute)
					else:
						date_minute = str(date_minute)
	
					# ========================================================
					# Read SEVIRI L1B and cloud type file, apply calibration coefficients, downscale CT
					# ========================================================

					sevl15filename = FindFile(l1b_inpath+date_year+"/"+date_month+"/"+date_day+"/", date_year+date_month+date_day+"t"+date_hour+date_minute, ".h5")
					sevctfilename = FindFile(ct_inpath, "S_NWC_CT_MSG2_rss-eu-VISIR_"+date_year+date_month+date_day+"T"+date_hour+date_minute, ".nc")
	
					if sevl15filename != "-1" and sevctfilename != "-1" and date_month != "08":
						# Read L1B
						Satellite_azimuth, Satellite_zenith, Sun_azimuth, Sun_zenith, Cal_Slope, Cal_Offset, Lambda_C, F_0, Ref_Slope, Ref_Offset, Nu_C, Alpha, Beta, imagech1, imagech2, imagech3, imagech4, imagech5, imagech6, imagech7, imagech8, imagech9, imagech10, imagech11, imagech12, PC_SouthLine_NONHRV, PC_NorthLine_NONHRV, PC_EastColumn_NONHRV, PC_WestColumn_NONHRV, PC_SouthLine_HRV, PC_NorthLine_HRV, PC_EastColumn_HRV, PC_WestColumn_HRV = SevReadL15fileTropos(sevl15filename)
		
						if PC_SouthLine_NONHRV != "-1":
		
							# Read CT
							sev_ct, sev_ctcumuli, sev_ctmulti, sev_ctconditions, sev_ctquality = SevReadCTfileTropos(sevctfilename)
							# Get CT flag values
							sev_ctconditionsflag, sev_ctqualityflag = SevCTFlags(sev_ctconditions, sev_ctquality)
		
							# Apply reflectance calibration
							sev_ref, sev_refhr = SevReflectancesTropos(imagech1, imagech2, imagech3, imagech12, Ref_Slope, Ref_Offset, Sun_zenith)
	
							# Reduce dataset to Germany domain if wanted
							if area == "DE":
								y0 = 60
								y1 = 300
								x0 = 320
								x1 = 720
			
								sy = slice(y0*3,y1*3)
								sx = slice(x0*3,x1*3)
								sev_refhr = np.array(sev_refhr[sy,sx])
								sy = slice(y0,y1)
								sx = slice(x0,x1)
								sev_ref = np.array(sev_ref[sy,sx,:])
								sev_ct = np.array(sev_ct[sy,sx])
								sev_ctconditionsflag = np.array(sev_ctconditionsflag[sy,sx,:])
			
							# Downscaling of sev_ct
							sev_ct_HRV = np.zeros((len(sev_refhr[:,0]),len(sev_refhr[0,:])), dtype=np.int)
							sev_ctconditionsflag_HRV = np.zeros((len(sev_refhr[:,0]),len(sev_refhr[0,:]),22), dtype=np.int)
							for i in range(0, len(sev_ref[:,0,0])):
								for j in range(0, len(sev_ref[0,:,0])):
									sev_ct_HRV[i*3:(i+1)*3,j*3:(j+1)*3] = sev_ct[i,j]
									sev_ctconditionsflag_HRV[i*3:(i+1)*3,j*3:(j+1)*3,:] = sev_ctconditionsflag[i,j,:]
		
							original_3km_cloudmask = np.zeros((len(sev_refhr[:,0]),len(sev_refhr[0,:])), dtype=np.int)
							original_3km_cloudmask[sev_ct_HRV >= 5] = 1 # we flag everything that is cloudy
			
							ffff = fff # this is a fake variable that we need for line 276-278
							tttt = tt_index # this is a fake variable that we need for line 276-278








							# ========================================================
							# Write the HRV reflectance, cloud type and cloud type conditions information into the daily arrays
							# ========================================================
							print fff,tt_index,np.float(date_hour)+np.float(date_minute)/np.float(60),csctime_array[ttt,0]+30/np.float(60),date_minute,tt
							sev_dailyhrv[:,:,fff,tt_index] = sev_refhr
							sev_dailyct[:,:,fff,tt_index] = sev_ct_HRV
							dummy = np.zeros((dimy,dimx), dtype=np.int)
							dummy[sev_ctconditionsflag_HRV[:,:,5]==1] = 1
							dummy[sev_ctconditionsflag_HRV[:,:,6]==1] = 2
							dummy[sev_ctconditionsflag_HRV[:,:,7]==1] = 3
							sev_dailyctc[:,:,fff,tt_index] = dummy
							sev_dailyoricm[:,:,fff,tt_index] = original_3km_cloudmask








				if np.float(date_hour)+np.float(date_minute)/np.float(60) == csctime_array[ttt,0]+30/np.float(60) or np.float(date_hour) == 16:
					# ========================================================
					# Determine thresholds
					# ========================================================
	
					# Define cloud type conditions variable with same dimensions as csc_data
					
					new_ctc_csc = np.zeros((dimy,dimx, len(csctime_array)), dtype=np.float32)
					for l in range(0, len(csctime_array)):
						dummy[:,:] = 0
						dummy[sev_dailyctc[:,:,ffff,tttt]==1] = 1
						dummy[sev_dailyctc[:,:,ffff,tttt]==2] = 2
						dummy[sev_dailyctc[:,:,ffff,tttt]==3] = 3
						new_ctc_csc[:,:,l] = dummy
	
					# Spatial and temporal median of csc_data
	
					landmedian_t_s = np.percentile(csc_data[(csc_data>0) & (new_ctc_csc==1)], 50)
					seamedian_t_s = np.percentile(csc_data[(csc_data>0) & (new_ctc_csc==2)], 50)
					coastmedian_t_s = np.percentile(csc_data[(csc_data>0) & (new_ctc_csc==3)], 50)
	
					# Temporal median of csc_data at each location
	
					landr_norm = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.float32)
					sear_norm = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.float32)
					coastr_norm = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.float32)
					landr_norm[:,:,:,:] = -999
					sear_norm[:,:,:,:] = -999
					coastr_norm[:,:,:,:] = -999
		
					for ii in range(0, len(csc_data[:,0])):
						for jj in range(0, len(csc_data[0,:])):
				
							median_t_dummy = np.zeros((len(csctime_array)), dtype=np.float32)
							median_t_dummy[:] = csc_data[ii,jj,:]
							median_t = 0
							if (median_t_dummy[median_t_dummy>0]).sum() > 0:
								median_t = np.percentile(median_t_dummy[median_t_dummy>0], 50)
						
							if  landmedian_t_s > 0 and median_t > 0 and new_ctc_csc[ii,jj,0] == 1:
								landr_norm[ii,jj,:,:] = csc_data[ii,jj,ttt] - (median_t - landmedian_t_s)
							if  seamedian_t_s > 0 and median_t > 0 and new_ctc_csc[ii,jj,0] == 2:
								sear_norm[ii,jj,:,:] = csc_data[ii,jj,ttt] - (median_t - seamedian_t_s)
							if  coastmedian_t_s > 0 and median_t > 0 and new_ctc_csc[ii,jj,0] == 3:
								coastr_norm[ii,jj,:,:] = csc_data[ii,jj,ttt] - (median_t - coastmedian_t_s)
	
					# Define variables to hold the HR cloud mask
	
					hr_mask = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.int)
					hr_mask_flags = np.zeros((dimy,dimx,julian_step, np.int(60./np.float(dt))), dtype=np.int)
					landthreshold[0,ttt] = np.nan
					seathreshold[0,ttt] = np.nan
					coastthreshold[0,ttt] = np.nan

					# Iterative determination of the HRV threshold via Matthews correlation coefficient for land surfaces
					if ((landr_norm>-999) & (sev_dailyctc==1) & (sev_dailyct != 3) & (sev_dailyoricm==1)).sum() > 0:
		
						# Range of the first guess of dummy threshold a
						dummy_a1 = 0.
						dummy_a2 = 0.5
					
						# Iterative calculation of best a, first in steps of 0.1, then in steps of 0.01 and then in steps of 0.001
						for n in range(1, 4):
							a_array = np.arange(dummy_a1,dummy_a2,1./(10.**n))
							min_val = -1
							for nn in range(0, len(a_array)):
								hr_mask[:,:,:,:] = 0
								hr_mask_flags[:,:,:,:] = 0
								
								hr_mask[(sev_dailyhrv > np.percentile(landr_norm[(landr_norm>-999) & (sev_dailyctc==1) & (sev_dailyct != 3) & (sev_dailyct<11)], 50) + a_array[nn]) & (landr_norm>-999) & (sev_dailyctc==1) & (sev_dailyct != 3) & (sev_dailyct<11)] = 1
								hr_mask_flags[(landr_norm>-999) & (sev_dailyctc==1) & (sev_dailyct != 3) & (sev_dailyct<11)] = 1
								
								matthews_corr, accuracy = MatthewsCorr(sev_dailyoricm, hr_mask, hr_mask_flags)
				
								if matthews_corr > min_val and accuracy > 0.2:
									dummy_a1 = a_array[nn]-1./(10.**n)
									dummy_a2 = a_array[nn]+1./(10.**n)
									min_val = matthews_corr
										
									landthreshold[0,ttt] = np.mean([dummy_a1,dummy_a2]) + np.percentile(landr_norm[(landr_norm>-999) & (sev_dailyctc==1) & (sev_dailyct != 3) & (sev_dailyct<11)], 50)

					# Iterative determination of the HRV threshold via Matthews correlation coefficient for sea surfaces
					if ((sear_norm>-999) & (sev_dailyctc==2) & (sev_dailyct != 4) & (sev_dailyoricm==1)).sum() > 0:

						# Range of the first guess of dummy threshold a
						dummy_a1 = 0.
						dummy_a2 = 0.5
		
						# Iterative calculation of best a, first in steps of 0.1, then in steps of 0.01 and then in steps of 0.001
						for n in range(1, 4):
							a_array = np.arange(dummy_a1,dummy_a2,1./(10.**n))
							min_val = -1
							for nn in range(0, len(a_array)):
								hr_mask[:,:,:,:] = 0
								hr_mask_flags[:,:,:,:] = 0
				
								hr_mask[(sev_dailyhrv > np.percentile(sear_norm[(sear_norm>-999) & (sev_dailyctc==2) & (sev_dailyct != 4) & (sev_dailyct<11)], 50) + a_array[nn]) & (sear_norm>-999) & (sev_dailyctc==2) & (sev_dailyct != 4) & (sev_dailyct<11)] = 1
								hr_mask_flags[(sear_norm>-999) & (sev_dailyctc==2) & (sev_dailyct != 4) & (sev_dailyct<11)] = 1
				
								matthews_corr, accuracy = MatthewsCorr(sev_dailyoricm, hr_mask, hr_mask_flags)
				
								if matthews_corr > min_val and accuracy > 0.2:
									dummy_a1 = a_array[nn]-1./(10.**n)
									dummy_a2 = a_array[nn]+1./(10.**n)
									min_val = matthews_corr
					
									seathreshold[0,ttt] = np.mean([dummy_a1,dummy_a2]) + np.percentile(sear_norm[(sear_norm>-999) & (sev_dailyctc==2) & (sev_dailyct != 4) & (sev_dailyct<11)], 50)

					# Iterative determination of the HRV threshold via Matthews correlation coefficient for coastal surfaces
					if ((coastr_norm>-999) & (sev_dailyctc==3) & (sev_dailyct != 3) & (sev_dailyct != 4) & (sev_dailyoricm==1)).sum() > 0:

						# Range of the first guess of dummy threshold a
						dummy_a1 = 0.
						dummy_a2 = 0.5
		
						# Iterative calculation of best a, first in steps of 0.1, then in steps of 0.01 and then in steps of 0.001
						for n in range(1, 4):
							a_array = np.arange(dummy_a1,dummy_a2,1./(10.**n))
							min_val = -1
							for nn in range(0, len(a_array)):
								hr_mask[:,:,:,:] = 0
								hr_mask_flags[:,:,:,:] = 0
				
								hr_mask[(sev_dailyhrv > np.percentile(coastr_norm[(coastr_norm>-999) & (sev_dailyctc==3) & (sev_dailyct != 3) & (sev_dailyct != 4) & (sev_dailyct<11)], 50) + a_array[nn]) & (coastr_norm>-999) & (sev_dailyctc==3) & (sev_dailyct != 3) & (sev_dailyct != 4) & (sev_dailyct<11)] = 1
								hr_mask_flags[(coastr_norm>-999) & (sev_dailyctc==3) & (sev_dailyct != 3) & (sev_dailyct != 4) & (sev_dailyct<11)] = 1
				
								matthews_corr, accuracy = MatthewsCorr(sev_dailyoricm, hr_mask, hr_mask_flags)
				
								if matthews_corr > min_val and accuracy > 0.2:
									dummy_a1 = a_array[nn]-1./(10.**n)
									dummy_a2 = a_array[nn]+1./(10.**n)
									min_val = matthews_corr
					
									coastthreshold[0,ttt] = np.mean([dummy_a1,dummy_a2]) + np.percentile(coastr_norm[(coastr_norm>-999) & (sev_dailyctc==3) & (sev_dailyct != 3) & (sev_dailyct != 4) & (sev_dailyct<11)], 50)
				
					# ========================================================
    				# Derive new cloud mask clear sky composite
    				# ========================================================
    				
    				# Step 1: set all reflectances to 0 where there are clouds, the proximity flag is 1
            		sev_dailyhrv[(sev_dailyctc == 1) & (sev_dailyhrv>landthreshold[0,ttt]) & (sev_dailyct<11)] = 0
            		sev_dailyhrv[(sev_dailyctc == 2) & (sev_dailyhrv>seathreshold[0,ttt]) & (sev_dailyct<11)] = 0
            		sev_dailyhrv[(sev_dailyctc == 3) & (sev_dailyhrv>coastthreshold[0,ttt]) & (sev_dailyct<11)] = 0
            		sev_dailyhrv[sev_dailyct>=11] = 0
            		
            		# Are there any HRV reflectances > 0 (at night unlikely)
            		if (sev_dailyhrv > 0).sum() > 0.:
            			# Step 2: derive the median of elements that have sev_dailyhrstatus = 1
            			for i in range(0, len(sev_refhr[:,0])):
            				for j in range(0, len(sev_refhr[0,:])):
            					dummy = np.zeros((julian_step,60./np.float(dt)), dtype=np.float32)
            					dummy[:,:] = sev_dailyhrv[i,j,:,:]
            					if (dummy > 0).sum() > 0:
            						csc_data[i,j,ttt] = np.percentile(dummy[dummy>0], 50)
				
					# ========================================================
    				# Calculate difference to old thresholds
    				# ========================================================
    				
    				it_loop_land = np.nanmean((landthreshold-landthreshold_old)**2)
    				it_loop_sea = np.nanmean((seathreshold-seathreshold_old)**2)
    				it_loop_coast = np.nanmean((coastthreshold-coastthreshold_old)**2)
    				
    				
    				
    				
    				
    				
    				
    				tt_index = -1
    				ttt += 1
    				sev_dailyhrv[:,:,:,:] = 0
    				sev_dailyct[:,:,:,:] = 0
    				sev_dailyctc[:,:,:,:] = 0
    				sev_dailyoricm[:,:,:,:] = 0
	
                tt += 1
	
                print "Finished calculating HR cloud mask thresholds for "+str(julian_step)+" day interval "+str(ff+1)+" of "+str(len(julian_array)-1)+", time step "+str(tt)+" of "+str(len(time_array))+" ("+str(date_hour)+":"+str(date_minute)+")"
                print int_loop,seathreshold








    # ========================================================
    # Write h5 file for each ff julian_array step, which includes the tt time_array steps
    # ========================================================
	if ff < len(julian_array)-2:
		sevclearfilename_dummy11, sevclearfilename_dummy12 = JulianDayToDate(date_year, julian_array[ff])
		sevclearfilename_dummy21, sevclearfilename_dummy22 = JulianDayToDate(date_year, julian_array[ff+1])
	else:
		sevclearfilename_dummy11, sevclearfilename_dummy12 = JulianDayToDate(date_year, julian_array[ff])
		sevclearfilename_dummy21, sevclearfilename_dummy22 = JulianDayToDate(date_year, julian_array[len(julian_array)-1])
	if sevclearfilename_dummy11 < 10:
		sevclearfilename_dummy11 = "0"+str(sevclearfilename_dummy11)
	else:
		sevclearfilename_dummy11 = str(sevclearfilename_dummy11)
	if sevclearfilename_dummy12 < 10:
		sevclearfilename_dummy12 = "0"+str(sevclearfilename_dummy12)
	else:
		sevclearfilename_dummy12 = str(sevclearfilename_dummy12)
	if sevclearfilename_dummy21 < 10:
		sevclearfilename_dummy21 = "0"+str(sevclearfilename_dummy21)
	else:
		sevclearfilename_dummy21 = str(sevclearfilename_dummy21)
	if sevclearfilename_dummy22 < 10:
		sevclearfilename_dummy22 = "0"+str(sevclearfilename_dummy22)
	else:
		sevclearfilename_dummy22 = str(sevclearfilename_dummy22)

	sevthreshfilename = thresh_outpath + "test/msg2-sevi-"+str(date_year)+"."+sevclearfilename_dummy11+sevclearfilename_dummy12+"-"+sevclearfilename_dummy21+sevclearfilename_dummy22+"-thresholds-rss-eu.c2.h5"
	SevWriteThreshFile(sevthreshfilename, landthreshold, seathreshold, coastthreshold, csctime_array, date_year, julian_step, dt)

print "Successfully derived HR cloud mask for "+str(date_year)
sys.exit("stop")

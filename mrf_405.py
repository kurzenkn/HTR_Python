# -*- coding: utf-8 -*-
"""
Created on Mon Jul 08 14:17:34 2019
Revised on Thu May 12 2022
******************************************************************************
language		: python
requirement		: a folder in the working directory called "hardwares" which contains the files: flames_s.py, arduino_control.py, multiplexer.py as indicated here:
				: https://github.com/LuWang04/python-mfr-feedback-control/tree/master/hardwares
authors			: Lu Wang, Chemical Engineering, University of Southern California
				: Ricki Chairil, Chemical Engineering, University of Southern California
function		: script performing routine IR and/or fluorescent spectral analysis and automatic stage movement of product photoluminescent properties
				  in CsPbBr3 millifluidic reactor
				: note pressure optimization will not be performed in this script. For pressure optimization, refer to mrf_feedback_control.py
******************************************************************************
"""

import datetime
import time
import threading
import logging
import numpy as np
import tkinter
from tkinter import filedialog
import os
import serial
import pandas as pd
from serial.tools.list_ports import comports
from pyfirmata import ArduinoMega 
from pyfirmata import util

# files created by LuWang
from hardwares import flame_s
from hardwares import arduino_control as ard_contr
from hardwares import mfcs
from hardwares import multiplexer as mux
# import nelder_mead as nm

# ====================== Check status of all the hardwares ===========================


# ------- locate the Arduino board for IR reading ------
ir_port = None
for comport in comports():
		SN = comport.serial_number
		if SN == '95530343235351605281':
			ir_port = comport.device # e.g. 'COM5'
if ir_port:
	print('The port of IR Arduino board (Mega_1) is :', ir_port)
else:
	input('The IR Arduino board (Mega_1) is not available. Press enter to continue!')


# ----- ensure linear stage is at initial position -----
stage_initial = input('Is linear stage at initial position? [y/n]')
if stage_initial in ['n','N','No','NO']:
	ard_contr.custom_move_stage()

# ------ initialize MFCS-EZ Fluigent pressure unit ------
#mfcs.init_mfcs() # 'MFCS is normal' means it's normal, START light is green on the unit
				 # 'MFCS is reset' means the valve is closed, STOP light is red
#read pressure in each channel
#mfcs.read_pressure()


# ============== choose a folder to save all the data ==================
root = tkinter.Tk()
root.withdraw() #use to hide tkinter window
file_dir = filedialog.askdirectory(parent=root, # initialdir=currdir, 
					title='Please select a directory to save data')
print(file_dir)


# ================ logging format ===================
logging.basicConfig(level=logging.DEBUG, 
	format='(%(threadName)-10s) %(message)s')


# ============== initialize optimization process ================
#method = nm.method
#x_list = nm.x_list
#f_list = nm.f_list
#reach_goal = False





# ################################## Function for reading IR signals #####################################
def ir_mux(ir_name):
	logging.debug('Starting')
	
	# ------- connect to IR hardwares --------
	board_port = ir_port # The port is from the USB port
	ir_board = ArduinoMega(board_port) # Creat the Serial port for ArduinoMega board
	# first multiplexer, define switch and signal pins on boards
	PinS, PinSig = mux.DefinePin(ir_board, [30,31,32,33],3)
	# second multiplexer, define switch and signal pins on boards
	PinS2, PinSig2 = mux.DefinePin(ir_board, [40,41,42,43],4)
	logging.debug('connected to IR board')

	# ------- initialize hardwares and variables ------
	# start switches with connectin to channel 0
	mux.SwitchMUX(PinS, 0)	
	mux.SwitchMUX(PinS2, 0)
	Data = {} # Creat an empty dictionary
	num_of_data = 2000

	# ------ start Iterator to avoid serial overflow -----
	it = util.Iterator(ir_board) # need this iteration otherwise pin just reporting 'None'
	it.start()
	time.sleep(0.5)

	# ------ start collect data ------
	for n in range(32):
		# logging.debug('IR channel '+str(n+1)+'...')
		print("IR channel "+str(n+1)+"...")
		t = []
		v = []
		t0 = time.time()
		if n<16:	# It's on the first mutiplexer
			mux.SwitchMUX(PinS,n)
			time.sleep(0.1)
			for m in range(num_of_data):
				t.append(time.time()-t0)
				v.append(float(mux.Read(PinSig))) # mux returns a string, need to transfer to float number
				time.sleep(0.000001)
			Data['Time'+str(n+1)] = t
			Data['Chan'+str(n+1)] = v
			
		else:		# switch to the second mutiplexer, need to re-define channel number
			m = n-16
			mux.SwitchMUX(PinS2,m)
			time.sleep(0.1)
			for m in range(num_of_data):
				t.append(time.time()-t0)
				v.append(float(mux.Read(PinSig2))) # MUX returns a string, need to transfer to float number
				time.sleep(0.000001)
			Data['Time'+str(n+1)] = t
			Data['Chan'+str(n+1)] = v
					
	ir_board.exit()
	np.save(os.path.join(file_dir,ir_name), Data) # use np.load(filename) to read data, and dict will be in .item()
	
	freq = []
	drops = []
	for i in range(32):
		time_title = 'Time'+str(i+1)
		chan_title = 'Chan'+str(i+1)
		time_stamp = Data[time_title]
		ir = Data[chan_title]
		ir = np.array(ir) # convert list to array, for the following process
		
		rel_ir = -1*ir+ max(ir)
		cut_height = 0.3* max(rel_ir)
		spec_cut_half = rel_ir -cut_height
		fwhm_index = np.where(np.diff(np.sign(spec_cut_half)))[0]
		peak_numbs = len(fwhm_index)/2
		# time lapse between reading is 1ms, the unit below is ms
		liq_total_time = sum(counter > cut_height for counter in rel_ir)
		
		freq.append(peak_numbs/(time_stamp[-1]-time_stamp[0])) # unit Hz
		drops.append(liq_total_time/peak_numbs) # unit: ms/drop
		
	print('########### Droplet frequency is  '+str(np.mean(freq))+'  Hz ##############')
	print('########### Droplet size is  '+str(np.mean(drops))+'  ms/drop ##############')

	logging.debug('Exiting')




# ###################################### PL-UV-vis spectra function #####################################
def pl_abs(pl_name, abs_name):
	logging.debug('Starting')
	
	peaks = [] 
	fwhms = [] 
	pl_info = {} 
	abs_info = {}

	# start to read spectra
	for chan in range(16):
							
		ard_contr.choose_chan(chan)
		
		print('=====================\nPL channel '+str(chan+1))
		ard_contr.light_source("LED")
		flame_s.get_PL(chan) # obtain PL spectral data for channel
		time.sleep(1)
		pl_time, pl_spec, peak_wave, fwhm = flame_s.read_pl()
		
		pl_info['Time'+str(chan+1)] = pl_time
		pl_info['Chan'+str(chan+1)] = pl_spec
		print('peak wave: ', peak_wave)
		print('FWHM: ', fwhm)

		ard_contr.light_source(None)
		time.sleep(1)

		print('Abs channel '+str(chan+1))
		ard_contr.light_source("UV-vis")
		flame_s.get_abs(chan) # obtain Abs spectral data for channel
		time.sleep(1)
		abs_time, abs_spec = flame_s.read_abs(chan)
		abs_info['Time'+str(chan+1)] = abs_time
		abs_info['Chan'+str(chan+1)] = abs_spec

		ard_contr.light_source(None)	# turn off both lights
		time.sleep(1)
		ard_contr.optical_switch(chan,0)	# let solenoid rest
		
		
		peaks.append(peak_wave) 
		fwhms.append(fwhm)
							
	
	np.save(os.path.join(file_dir,pl_name), pl_info) # load with np.load(filename), dict will be in .item()
	np.save(os.path.join(file_dir,abs_name), abs_info) # load with np.load(filename), dict will be in .item()
	np.savetxt(os.path.join(file_dir,peak_name), peaks)
	np.savetxt(os.path.join(file_dir,fwhm_name), fwhms)

	logging.debug('Exiting')






# ######################################### Main script #############################################
if __name__ == "__main__":
	try:
		# [p1_set,p3_set,p4_set] = nm.x_list[0]
		
		# mfcs.set_pressure(1,p1_set)		# P[Cs-Pb]
		# mfcs.set_pressure(3,p3_set)		# P[gas]
		# mfcs.set_pressure(4,p4_set)		# P[Br]
		
		
		time.sleep(60) # wait for liquid to reach detector modules
			
		for _ in range(1): # number of passes to make
			t_now = datetime.datetime.now()
			
			# assign file names
			time_str = str(t_now.hour).zfill(2)+str(t_now.minute).zfill(2)
			npy_extension = time_str+'_'+"PA=300"+'_'+"PB=200"+'_'+"PC=300"+'.npy' # name files according to MFCS pressures used; convert to .npy array 
			out_extension = time_str+'_'+"PA=300"+'_'+"PB=200"+'_'+"PC=300"+'.out'  # name files according to MFCS pressures used; convert to .out file
			pl_name = 'PL_'+ npy_extension
			abs_name = 'Abs_'+ npy_extension 
			ir_name = 'IR_'+ npy_extension
			peak_name = 'peaks_'+ out_extension
			fwhm_name = 'fwhms_1'+ out_extension
			
			pl_abs(pl_name, abs_name) # call this alone if IR not used
			'''
			# activate lines 254-264 if multithreading IR module as well
			# define then start ir_mux thread and pl_abs thread
			t1 = threading.Thread(name = 'IR', target = ir_mux, args = ([ir_name]))
			t2 = threading.Thread(name = 'PL-UV-vis', target = pl_abs, args = (pl_name, abs_name))
				
			t1.start()
			t2.start()
			
			t1.join()
			print(threading.enumerate())
			t2.join()
			print(threading.enumerate())
			'''
			# process PL data in real-time					
			fwhms = np.loadtxt(os.path.join(file_dir,fwhm_name))
						
			fwhm_nansize = (np.isnan(fwhms)).sum(0) # count numbers of fwhms that are 'nan'
			fwhm_size = len(fwhms) - fwhm_nansize
			if fwhm_nansize > 2:
				fwhms = [60 if np.isnan(x) else x for x in fwhms]
				fwhm_size = 16
			averg_width = np.nanmean(fwhms)
			stand_err = np.nanstd(fwhms)/np.sqrt(fwhm_size)
			confidence = averg_width + 1.96 * stand_err # upper endpoint of 95% confidence interval
			
			
			print('##################### FWHM is: '+ str(averg_width)+'+/-'+str(stand_err) +' ######################')
			print('confidence: ', confidence)
			'''
			if confidence >35: 
				if reach_goal == True: # had previously reached goal, should start to optimize again
					reach_goal = False
					if confidence > 37.5: # big deviation, probably confronting a disturbance, start a new optimation
						method = None
						# define a new nm.x_list based on current pressures
						nm.new_x_list(p1_set, p3_set, p4_set, confidence)
					# else: # resume the optimization without setting the new nm.x_list, as this is just a small deviation.
	
				method, [p1_set,p3_set,p4_set] = nm.simplex(method,confidence)
				mfcs.set_pressure(1, p1_set)
				mfcs.set_pressure(3, p3_set)
				mfcs.set_pressure(4, p4_set)
				
				# move stage back to position 0 and wait till pressures stablize.
				print('move back to position 0')
				ard_contr.stage_return()
				time.sleep(30)
			else:
				if min(p1_set, p4_set) < 398: # increase all the pressures to increase throughput
					[p1_set,p3_set,p4_set] = nm.expand_pressure(p1_set,p3_set,p4_set,400)
					mfcs.set_pressure(1, p1_set)
					mfcs.set_pressure(3, p3_set)
					mfcs.set_pressure(4, p4_set)
					# move stage back to position 0.
					print('move back to position 0')
					ard_contr.stage_return()
					reach_goal = True
					print('reach_goal? ', reach_goal)
					# fwhm narrow enough, can decrease the frequency of detection.
					print('## desired product, but need to increase pressure ##')
					time.sleep(30)
				else:
					# move stage back to position 0, must wait after 16 channel processes are done.
					print('move back to position 0')
					ard_contr.stage_return()
					reach_goal = True
					print('reach_goal? ', reach_goal)
					# fwhm narrow enough, can decrease the frequency of detection.
					print('######################### desired product, wait for 60s till next detection ##########################')
					time.sleep(60)
		'''	
				

	except KeyboardInterrupt:
		print('Program terminated by key interruption')
		
	finally:
		print('Finally module')

		# export intensities as csvs (this will only do one iteration)
		# spectrum = flames.intensities()
		# df = pd.DataFrame(spectrum)
		# df.to_csv("PL_spectra.csv") # print intensity data for blank to a csv file via Pandas
		
		# ard_contr.exit_boards()
		# mfcs.exit_mfcs()
		# ir_board.exit()
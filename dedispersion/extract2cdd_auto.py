import numpy as np
import matplotlib.pyplot as plt
import fnmatch
import sys
import os
import pandas as pd
import csv
import itertools
sys.path.insert(0, os.path.abspath('../extractor'))

def _read_data(database):
	
	'''Reads in database (.csv format) containing file 
	paths for fil and raw files'''

	filepaths = pd.read_csv(str(database))
	filpaths = filepaths.iloc[:,0]
	rawpaths = filepaths.iloc[1, :][1:]
	fieldnames = filepaths.columns[1:]
	
	#If more than one csv is being analyzed, try this:
		#csvs = [str(csv)]
		#for csv in os.listdir(str(csv_dir)):
		#	if csv.endswith('.csv'):
		#		csvs.append(csv)

	return filepaths, filpaths, rawpaths, fieldnames


def _parse_spandak(csvs, sd_grade, DM_min=100, DM_max=2000, intervene=False):

	'''Filters each file for desired DM and SNR limits, and finally duplicate 
	detections by SPANDAK. Elements offered by the SPANDAK-generated
	csv are then parsed to extract key information about each burst 
	having a grade B or (not and) C'''

	for csv in csvs:
		
		#STAGE 1: FILTER OUT CANDIDATES THAT FALL OUTSIDE DESIRED DM RANGE

		#Read in .csv produced by SPANDAK
		cands = pd.read_csv(csv)

		#Isolate SPANDAK candidates (probable if *not* using ML functionality)
		sd_idx_nosnrfil_nodmfil = cands[cands.iloc[:, :]['Category']==str(sd_grade)].index.values

		#Exit if no SPANDAK candidates were detected
		if intervene == False:
			if len(sd_idx_nosnrfil_nodmfil) == 0:
				sys.exit('Sorry No ' + str(sd_grade) + ' Candidates Were Detected, Darn!')		

		#Find DMs for SPANDAK candidates
		DMs_nosnrfil_nodmfil = [[i for i in cands.loc[:, :]['DM']][b] for b in \
					sd_idx_nosnrfil_nodmfil]

		#Find DM-thresholded SPANDAK candidates
		sd_idx_nosnrfil = [sd_idx_nosnrfil_nodmfil[i] for i in \
			np.arange(len(sd_idx_nosnrfil_nodmfil)) if \
			float(DM_min) < DMs_nosnrfil_nodmfil[i] < float(DM_max)]

		#Exit if no bursts fall within desired dm range
		if intervene == False:	
			if len(sd_idx_nosnrfil) == 0:
				sys.exit(str(sd_grade) + ' Candidates Were Detected, But Not Within The Right DM Range!')			

		#Find DMs for DM-thresholded SPANDAK candidates
		DMs_nosnrfil = [[i for i in cands.loc[:, :]['DM']][b] for b in \
					sd_idx_nosnrfil]

		#Find time centers for DM-thresholded SPANDAK candidates
		parse_center_nosnrfil = [[i for i in [j.split('_') for j in \
		cands.iloc[:, :]['PNGFILE']]][b][2] for b in sd_idx_nosnrfil]
		filenumber = [i for i in cands.loc[:, :]['filename']][0][-27:-25]

		#Correct for observation time offset in TOAs
		toas_nosnrfil = [(np.float(m[:-3]) + ((int(filenumber) - 11) * 1800)) for m in \
			parse_center_nosnrfil]

		#STAGE 2: IDENTIFY DUPLICATES BASED ON TOA PROXIMITY
		
		class Delta:
			def __init__(self, delta):
				self.last = None
				self.delta = delta
				self.key = 1
			def __call__(self, value):
				if self.last is not None and abs(self.last - value) > self.delta:
					# Compare with the last value (`self.last`)
					# If difference is larger than 0.02, advance to next project
					self.key += 1
				self.last = value  # Remeber the last value.
				return self.key

		#For candidates separated by less than 0.02 seconds (TOAs)
		dup_times = []
		for key, grp in itertools.groupby(np.sort(toas_nosnrfil), key=Delta(0.02)):
			dup_times.append(list(grp))

		#STAGE 3: SAVE MAX SNR CANDIDATE FOR IDENTIFIED DUPLICATES AND FILTER OUT THE REST

		#Find SNRs for all candidates
		SNRs = [i for i in cands.loc[:, :]['SNR']]

		#Fnd max SNR candidates and remove the rest
		for sublist in dup_times:
			if len(sublist) > 1:
				temp_SNRs = []
				for i in sublist:
					temp_SNRs.append(SNRs[np.where(toas_nosnrfil == i)[0][0]])
				max_snr = sublist[np.where(temp_SNRs == np.max(temp_SNRs))[0][0]]
				dup_times.append([max_snr])
				dup_times.remove(sublist)

		#Locate candidates which are fully filtered and unique
		sd_idx = []
		for dup in dup_times:
			sd_idx.append(sd_idx_nosnrfil[np.where(dup[0] == toas_nosnrfil)[0][0]])

		#Identify parameters for fully filtered and unique candidates

		parse_center = [[i for i in [j.split('_') for j in \
		cands.iloc[:, :]['PNGFILE']]][b][2] for b in sd_idx]
		filenumber = [i for i in cands.loc[:, :]['filename']][0][-27:-25]
		
		toas = [(np.float(m[:-3]) + ((int(filenumber) - 11) * 1800)) for m \
				in parse_center]

		#Find DMs for fully filtered unique candidates
		DMs = [[i for i in cands.loc[:, :]['DM']][b] for b in \
				sd_idx]

		#Find filenames for fully filtered unique candidates
		files = [[i for i in cands.loc[:, :]['filename']][b] for b in \
				sd_idx]

		#Find source names for fully filtered unique candidates
		sourcename = [[i for i in cands.loc[:, :]['SourceName']][b] \
						for b in sd_idx]
		
		#Time widths for fully filtered unique candidates
		time_widths = [[i for i in cands.loc[:, :]['WIDTH']][b] for b in \
						sd_idx]

		#STAGE 4: MANUAL INTERVENTION IN CASE BURSTS ARE TOO WEAK

		if intervene == True:
			toas = [18.4463]
			DMs = [547.3]
			start_times = [18.2463]
			end_times = [18.8463]

			'''
			toas = [18.4463, 26.4042, 224.2854, 263.3902, 263.4046, 277.3659, \
			315.0407, 323.3557, 344.7733, 560.8777, 580.6574, 597.6409, 652.5922, \
			26.4042, 691.0704, 704.1022, 26.4042, 963.5013, 993.3134, 1036.4929, \
			1071.8031, 1073.1521, 263.3902, 1440.9448, 1454.6557, 1613.6106, 1630.9135, \
			1896.9579, 1923.0727, 2170.2814, 263.4046, 2485.9404, 3311.8506, 3520.9209, \
			26.4042, 597.6409, 1454.6557, 560.8777, 1071.8031, 704.1022, 8933.081, 9462.0387, \
			9462.0414, 12932.3224, 14232.5532, 1036.4929, 15620.6784, 15758.9369, 16074.1381]
			DMs = [547.3, 646.1, 607.0, 592.2, 631.4, 636.6, 589.3, 585.3, 590.7, 581.2, 557.9, \
			574.4, 578.0, 646.1, 559.7, 575.9, 646.1, 566.9, 566.2, 565.7, 595.8, 596.0, 592.2, \
			622.9, 563.6, 618.4, 569.7, 584.0, 572.0, 646.8, 631.4, 569.0, 593.5, 581.3, 646.1, \
			574.4, 563.6, 581.2, 595.8, 575.9, 581.4, 582.1, 576.2, 609.0, 622.9, 565.7, 595.5, 579.2, 556.8]
			start_times = [18.2463, 26.2042, 224.08540000000002, 263.1902, 263.2046, 277.1659, \
			314.8407, 323.1557, 344.5733, 560.6777, 580.4574, 597.4408999999999, 652.3922, 26.2042, \
			690.8703999999999, 703.9022, 26.2042, 963.3013, 993.1134, 1036.2929, 1071.6031, 1072.9521, \
			263.1902, 1440.7448, 1454.4557, 1613.4106, 1630.7135, 1896.7579, 1922.8727, 2170.0814, \
			263.2046, 2485.7404, 3311.6506000000004, 3520.7209000000003, 26.2042, 597.4408999999999, \
			1454.4557, 560.6777, 1071.6031, 703.9022, 8932.881, 9461.838699999998, 9461.8414, \
			12932.122399999998, 14232.3532, 1036.2929, 15620.4784, 15758.7369, 16073.9381]	
			end_times = [18.8463, 26.804199999999998, 224.68540000000002, 263.79019999999997, \
			263.8046, 277.7659, 315.4407, 323.7557, 345.1733, 561.2777, 581.0574, 598.0409, \
			652.9922, 26.804199999999998, 691.4703999999999, 704.5022, 26.804199999999998, \
			963.9013, 993.7134, 1036.8929, 1072.2031000000002, 1073.5521, 263.79019999999997, \
			1441.3448, 1455.0557000000001, 1614.0106, 1631.3135000000002, 1897.3579000000002, \
			1923.4727, 2170.6814, 263.8046, 2486.3404, 3312.2506000000003, 3521.3209, 26.804199999999998, \
			598.0409, 1455.0557000000001, 561.2777, 1072.2031000000002, 704.5022, 8933.481, 9462.438699999999, \
			9462.4414, 12932.722399999999, 14232.9532, 1036.8929, 15621.0784, 15759.3369, 16074.5381]
			'''

		#Calculate dispersion delay in s for fully filtered unique candidates
		band = [[i for i in cands.loc[:, :]['BANDWIDTH']][b] for b in \
				sd_idx]
		tau_disp = []
		for B in np.arange(len(sd_idx)):	
			tau_disp.append((4.149e+3) * DMs[B] * ((band[B])**(-2)))

		#Calculate start and end times for raw voltage extraction

		#Time limits: widths adjusted to 3 * dispersion delay (i.e., tau_disp)
		if intervene == False:	
			start_times = [toas[b]-(1*tau_disp[b]) for b in \
					np.arange(len(sd_idx))]
		if intervene == False:	
			end_times = [toas[b]+(2*tau_disp[b]) for b in \
					np.arange(len(sd_idx))]


		#CHECK
		#Print candidate indices
		print(str(sd_grade) + ' TOAs: ', toas)
		print(str(sd_grade) + ' DMs: ', DMs)
		print(str(sd_grade) + ' Start Times: ', start_times)
		print(str(sd_grade) + ' End Times: ', end_times)

		return sd_idx, files, DMs, sourcename, time_widths, toas, \
				start_times, end_times, tau_disp, csv


def _extract_auto(rawpaths, fieldnames, sd_idx, files, filepaths, \
				start_times, end_times, csv, mkdir=False, extract=False):

	'''Constructs extraction command for raw voltages and stores information for
	plotting. Sections for subbanded extraction are commented out, but still 
	categorize candidates. This currently only handles 4 subbands'''

	#Store extraction commands in list
	extract_run_commands = []

	#For subbanded *extraction*, not categorization

	'''
	extract_run_commands = {}
	extract_run_commands['0'] = []
	extract_run_commands['1'] = []
	extract_run_commands['2'] = []
	extract_run_commands['3'] = []
	extract_run_commands['all'] = []
	'''

	'''Group together candidates with start times within a standard... 
	deviation of 0.2...works for now, can change..'''

	st_sorted = np.sort(start_times)
	gaps = [y - x for x, y in zip(st_sorted[:-1], st_sorted[1:])]
	sd = np.std(gaps)
	group_starts = [[[start_times.index(st_sorted[0]), st_sorted[0]]]]
	for t in st_sorted[1:]:
		if (t - group_starts[-1][-1][-1]) > 0.2:
			group_starts.append([])
		group_starts[-1].append([start_times.index(t), t])

	'''Store candidate idxs into subcategories within dictionary based
	on the subbands in which they were flagged by SPANDAK'''

	sub_cands = {}
	sub_cands['cand_0'] = []
	sub_cands['cand_1'] = []
	sub_cands['cand_2'] = []
	sub_cands['cand_3'] = []
	sub_cands['all'] = []
	sub_cands['combined'] = []
       
	#Store indices for subbanded candidates
	for B in np.arange(len(sd_idx)):
		if 'diced_0' in files[B]:
			sub_cands['cand_0'].append(B)
			sub_cands['combined'].append(B)
		elif 'diced_1' in files[B]:
			sub_cands['cand_1'].append(B)
			sub_cands['combined'].append(B)
		elif 'diced_2' in files[B]:
			sub_cands['cand_2'].append(B)
			sub_cands['combined'].append(B)
		elif 'diced_3' in files[B]:
			sub_cands['cand_3'].append(B)
			sub_cands['combined'].append(B)
		else:
			sub_cands['all'].append(B)
			sub_cands['combined'].append(B)

	#Store band information for plotting
	plot_bands = {}
	plot_bands['3.8_5.1'] = []
	plot_bands['3.8_6.4'] = []
	plot_bands['3.8_7.7'] = []
	plot_bands['5.3_6.4'] = []
	plot_bands['5.3_7.7'] = []
	plot_bands['5.3_9'] = []
	plot_bands['6.6_7.7'] = []
	plot_bands['6.6_9'] = []
	plot_bands['7.9_9'] = []
	plot_bands['3.8_9'] = []


	'''Determine band information by checking to see how many subbands a candidate has
	been found in (by TOA proximity). If a candidate is found in two adjacent subbands
	then the complete range of both subbands together is stored for plotting.'''
	for cgroup in group_starts:
		if len(cgroup) == 3:
			if cgroup[0][0] in sub_cands['all']:
				plot_bands['3.8_9'].append([cgroup[0][0], \
					start_times[cgroup[0][0]], end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] \
				in sub_cands['cand_2'] and \
					cgroup[2][0] in sub_cands['cand_1']:
				plot_bands['3.8_7.7'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
					sub_cands['cand_2'] and \
					cgroup[2][0] in sub_cands['cand_0']:
				plot_bands['3.8_9'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
					sub_cands['cand_1'] and \
					cgroup[2][0] in sub_cands['cand_0']:
				plot_bands['3.8_9'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2'] and cgroup[1][0] in \
					sub_cands['cand_1'] and \
					cgroup[2][0] in sub_cands['cand_0']:
				plot_bands['5.3_9'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
		if len(cgroup) == 2:
			if cgroup[0][0] in sub_cands['all']:
				plot_bands['3.8_9'].append([cgroup[0][0], \
					start_times[cgroup[0][0]], end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
				sub_cands['cand_2']:
				plot_bands['3.8_6.4'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
				sub_cands['cand_1']:
				plot_bands['3.8_7.7'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
				sub_cands['cand_0']:
				plot_bands['3.8_9'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2'] and cgroup[1][0] in \
				sub_cands['cand_1']:
				plot_bands['5.3_7.7'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2'] and cgroup[1][0] in \
				sub_cands['cand_0']:
				plot_bands['5.3_9'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_1'] and cgroup[1][0] in \
				sub_cands['cand_0']:
				plot_bands['6.6_9'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
		else:
			if cgroup[0][0] in sub_cands['all']:
				plot_bands['3.8_9'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3']:
				plot_bands['3.8_5.1'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2']:
				plot_bands['5.3_6.4'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_1']:
				plot_bands['6.6_7.7'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_0']:
				plot_bands['7.9_9'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
   
	#IF EXTRACTING SUBBAND VOLTAGES: CONSTRUCT EXTRACTION COMMANDS

	'''			
	for B in sub_cands['cand_0']:
		for raw0 in np.arange(0, 7):
			extract_run_0 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
			+ rawpaths[raw0] + ' ' + 'blc' + str(fieldnames[raw0][3:]) + files[B][75:-25] + ' ' \
			+ str(start_times[B]) + ' ' + str(end_times[B]) \
			+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
			+ str(start_times[B]) + '_' + str(end_times[B]) + '_7.9_9/'
			extract_run_commands['0'].append(extract_run_0)
	for B in sub_cands['cand_1']:
		for raw1 in np.arange(7, 14):
			extract_run_1 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
			+ rawpaths[raw1] + ' ' + 'blc' + str(fieldnames[raw1][3:]) + files[B][75:-25] + ' ' \
			+ str(start_times[B]) + ' ' + str(end_times[B]) \
			+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
			+ str(start_times[B]) + '_' + str(end_times[B]) + '_6.6_7.7/'
			extract_run_commands['1'].append(extract_run_1)
	for B in sub_cands['cand_2']:
		for raw2 in np.arange(14, 21):
			extract_run_2 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
			+ rawpaths[raw2] + ' ' + 'blc' + str(fieldnames[raw2][3:]) + files[B][75:-25] + ' ' \
			+ str(start_times[B]) + ' ' + str(end_times[B]) \
			+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
			+ str(start_times[B]) + '_' + str(end_times[B]) + '_5.3_6.4/'
			extract_run_commands['2'].append(extract_run_2)
	for B in sub_cands['cand_3']:
		for raw3 in np.arange(21, 29):
			extract_run_3 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
			+ rawpaths[raw3] + ' ' + 'blc' + str(fieldnames[raw3][3:]) + files[B][75:-25] + ' ' \
			+ str(start_times[B]) + ' ' + str(end_times[B]) \
			+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
			+ str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_5.1/'
			extract_run_commands['3'].append(extract_run_3)
	for B in sub_cands['combined']:
		for raw in np.arange(len(rawpaths)):
			extract_run = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
			+ rawpaths[raw] + ' ' + 'blc' + str(fieldnames[raw][3:]) + files[B][33:-25] + ' '  \
			+ str(start_times[B]) + ' ' + str(end_times[B]) \
			+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' +  str(start_times[B]) + '_' \
			+ str(end_times[B]) + '_3.8_9/raws/'
			extract_run_commands.append(extract_run)
	'''

	#Construct extraction command (this has been hard coded for R3 and 121102)
	#for B in sub_cands['all']:
	for B in np.arange(len(start_times)):										   #FOR MANUAL INTERVENTION
		for raw in np.arange(len(rawpaths)):
			extract_run = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
			+ rawpaths[raw] + ' ' + 'blc' + str(fieldnames[raw][3:]) + files[0][85:-25] + ' '  \
			+ str(start_times[B]) + ' ' + str(end_times[B]) \
			+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/'  \
			+  str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9/raws/'
			extract_run_commands.append(extract_run)
			#Substitute in when necessary
			#+ rawpaths[raw] + ' ' + 'blc7' + str(fieldnames[raw][4:]) + filepaths.iloc[:,0][1][5:-22] + ' ' \
			#+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/0003_raws/' \

	#IF EXTRACTING SUBBAND VOLTAGES: MAKES DIRECTORIES FOR EACH SUBBAND FOR STORING VOLTAGES

	'''
	sb_7_9_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) \
				+ '_' + str(end_times[B]) + '_7.9_9'
	sb_6_7_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) \
				+ '_' + str(end_times[B]) + '_6.6_7.7'
	sb_5_6_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) \
				+ '_' + str(end_times[B]) + '_5.3_6.4'
	sb_3_5_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) \
				+ '_' + str(end_times[B]) + '_3.8_5.1'
	
	for B in sub_cands['all']:
	#for B in np.arange(len(start_times)):	#FOR MANUAL INTERVENTION
	
		if not os.path.exists(sb_7_9_dir):
			os.mkdir(sb_7_9_dir)
		if not os.path.exists(sb_6_7_dir):
			os.mkdir(sb_5_6_dir)
		if not os.path.exists(sb_5_6_dir):
			os.mkdir(sb_5_6_dir)
		if not os.path.exists(sb_3_5_dir):
			os.mkdir(sb_3_5_dir)
	'''

	#For R3

	'''
	full_band_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' \
				+ str(csv[13:17]) + '/' + str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9'

	for B in sub_cands['all']:	
		if not os.path.exists(full_band_dir):
			os.mkdir(full_band_dir)
	'''

	#Form directories for storing raw voltages
	if mkdir == True:
		#Construct time directories
		time_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' \
					+ str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9'

		#Construct raw voltage directories
		raws_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' \
					+ str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9/raws'

		#for B in sub_cands['all']:
		for B in np.arange(len(start_times)): #FOR MANUAL INTERVENTION
			if not os.path.exists(time_dir):
				os.mkdir(time_dir)
			if not os.path.exists(raws_dir):
				os.mkdir(raws_dir)

	#Run extraction commands
	if extract == True:
		for erc in extract_run_commands:
			#print("Extraction Commands: ", erc)
			os.system(erc)

	return extract_run_commands, sub_cands, plot_bands


def _splice_auto(sub_cands, files, start_times, end_times, \
	ex_raws_path='/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/', splice=False):

	'''Splices extracted raw voltages'''

	#Parse and Form Raw File Splicing Commands
	splicer_run_commands = []

	#for B in sub_cands['combined']:
	for B in np.arange(len(start_times)):	#FOR MANUAL INTERVENTION	
		splicer_run = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/splicer_raw.py ' \
		+  str(ex_raws_path) + str(start_times[B]) + '_' + str(end_times[B]) \
		+ '_3.8_9/raws ' + '2 ' + 'spliced_' +  str(start_times[B]) + '_' + str(end_times[B]) + '.raw'
		splicer_run_commands.append(splicer_run)

	#Splice Extracted Raw Files Into Contiguous Raw File
	if splice == True:
		for src in splicer_run_commands:
			#print('Splice Raw Commands :', src)
			os.system(src)

	return splicer_run_commands

def _gen_par(sourcename, sd_idx, DMs, write_par=False):

	'''Writes par files'''

	source_int = [str(sub.split('_')[1])[3:] for sub in sourcename]

	#FRB121102 rotation frequency
	F0 = 47.680639719963075
	
	#R3 rotation frequency
	#F0 = 12

	par_file = []
	#for B in np.arange(len(sd_idx)):
	for B in np.arange(len(start_times)):  #FOR MANUAL INTERVENTION
		par_txt = 'PSR  ' + str(121102) + '\n' \
		+ 'RAJ  ' + '05:31:58.70' '\n' + 'DECJ  ' + '+33:08:52.5' + '\n' \
		+ 'C The following is (1+v/c)*f_topo where f_topo is 1/(2048*dt)' + '\n' \
		+ 'F0  ' + str(F0) + '\n' + 'DM  ' + str(DMs[B]) + '\n' + 'PEPOCH  ' \
		+ '57962.373622685185186' + '\n' + 'CLOCK  ' + 'UNCORR'
		par_file.append(par_txt)
		
	par_fil_paths = []
	#for B in np.arange(len(sd_idx)):
	for B in np.arange(len(start_times)): #FOR MANUAL INTERVENTION
		par_fil = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/pars/' \
				+ 'FRB_' + str(toas[B]).split('.')[0] + '.par'
		#par_fil = '/Users/jakobfaber/Documents/spandak_extended/SPANDAK_extension/pipeline_playground/parfiles/' \
				#+ 'FRB_' + str(source_int[B]) + '_' + str(sd_idx[B]) + '.par'
		par_fil_paths.append(par_fil)
	
	#Write par files
	if write_par == True:
		for B in np.arange(len(sd_idx)):
		#for B in np.arange(len(start_times)):	
			par = open(par_fil, "w")
			par.write(par_file[B])
			par.close()

	return source_int, par_file, par_fil_paths


def _cdd_auto(sub_cands, files, par_fil_paths, start_times, end_times, cdd=False):

	'''Performs coherent dedispersion with DSPSR on spliced raw voltages'''

	#Specify CDD parameters
	cepoch = 58178
	polar = 4
	phasebin = 2048
	p = 0
	chan = 9728
	samples = 1024 #number of MB

	#Parse and Form Coherent Dispersion Commands
	cdd_run_commands = []

	for B in np.arange(len(start_times)):
		cdd_run = 'dspsr ' + '-U ' + str(samples) + ' -F ' + str(chan) + ':D ' \
		+ ' -K ' + ' -d ' + str(polar) + ' -b  ' + str(phasebin) + ' -E ' \
		+ par_fil_paths[B] + ' -s -a psrfits -e fits ' \
		+ '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' \
		+ str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9/raws/' + 'spliced_' \
		+ str(start_times[B]) + '_' + str(end_times[B]) + '.raw'
		cdd_run_commands.append(cdd_run)

	#Coherently Dedisperse Spliced Raw Voltages
	if cdd == True:

		for cdd in cdd_run_commands:
			
			#Make directies for raw voltages (if necessary -- probably won't need it), but also for fits file storage
			time_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' + str(cdd.split('/')[16])
			fits_dir = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' + str(cdd.split('/')[16]) + '/fits'
			#print('Time', time_dir)
			#print('Fits', fits_dir)

			if not os.path.exists(time_dir):
				os.mkdir(time_dir)
			if not os.path.exists(fits_dir):
				os.mkdir(fits_dir) 

			#Funnel dspsr output into fits file directory made previously
			os.chdir(fits_dir)
			print('DSPSR Output Funnelling Into: ' + os.getcwd())
			os.system(cdd)
			print('Coherent Dedispersion Complete')



	return cdd_run_commands

if __name__ == "__main__":

	#Define input arguments: database, SPANDAK csv, SPANDAK grade
	database = sys.argv[1]
	csvs = [sys.argv[2]]
	sd_grade = sys.argv[3]
	
	#Read Database
	filepaths, filpaths, rawpaths, fieldnames = _read_data(database)
	#Identify Relevant Candidate Parameters
	sd_idx, files, DMs, sourcename, time_widths, toas, start_times, end_times, \
		tau_disp, csv = _parse_spandak(csvs, sd_grade, DM_min=320, DM_max=380, intervene=True)
	#Extract Raw Voltages
	extract_run_commands, sub_cands, plot_bands = _extract_auto(rawpaths, \
		fieldnames, sd_idx, files, filepaths, start_times, end_times, csv, mkdir=True, extract=True)
	#Splice Raw Voltages
	splicer_run_commands = _splice_auto(sub_cands, files, start_times, end_times, splice=True)
	#Generate Par Files
	source_int, par_file,par_fil_paths = _gen_par(sourcename, sd_idx, DMs, write_par=True)
	#Coherently Dedisperse Spliced Raw Voltages
	cdd_run_commands = _cdd_auto(sub_cands, files, par_fil_paths, start_times, end_times, cdd=True)









import datetime

import calculations as calc

filepath = "../measurements/"

columnnames_roi = {"BEGIN (ADC CH.)":"Begin", "END (ADC CH.)":"End", "UNCERTAINITY":"Uncertainty", "CENTROID (ADC CH.)":"Centroid", "FWHM (ADC CH.)":"Fwhm", "RANGE (ADC CH.)":"Range", "GROSS":"Gross", "NET":"Net", "PEAK":"Peak"}

#Filename of measurement for
#Source: "Co, Ba, Bi, Rn, Background"
#Region of Interest set? True/False
#Duration of measurement: 2, 4 (hours)
def filename(source, roi=False, duration=2, custom_postfix=""):
	"""Filename of measurement for
Source: "Co, Ba, Bi, Rn, Background"
Region of Interest set? True/False
Duration of measurement: 2, 4 (hours)"""
	fname = filepath
	calibration = False
	background = False

	if (source == "Co"):
		fname += "60_Co"
		calibration = True
	elif (source == "Ba"):
		fname += "133_Ba"
		calibration = True
	elif (source == "Bi"):
		fname += "207_Bi"
		calibration = True
	elif (source == "Rn"):
		fname += "Rn"
	elif (source == "Background"):
		fname += "UntergrundJuni2017"
		background = True

	if (calibration):
		fname += "_Calibration"
	elif (not background):
		if (duration == 2):
			fname += "_2h"
		elif (duration == 4):
			fname += "_4h"
		else:
			fname += "_invalid_time"
		
		fname += "_measurement"
	
	if (roi):
		fname += "_ROI"

	fname += custom_postfix
	fname += ".txt"

	return fname


#Cuts line in words separated by more than 2 spaces removing those
def parse_line(line, numbers=False):
	"""Cuts line in words separated by more than 2 spaces removing those"""
	wordlist = []
	numberlist = []
	returnlist = []

	while (line.find("   ") != -1):	
		line = line.replace("   ", "  ")
	wordlist = line.split("  ")

	#Remove "\r\n" from end of last item
	wordlist[len(wordlist)-1] = wordlist[len(wordlist)-1].split("\r")[0]
	#Remove first item ("")
	wordlist.pop(0)

	if (numbers):
		for word in wordlist:
			try:
				numberlist.append((float)(word))
			except ValueError:
				tmplist = []
				for part in word.split("-"):
					tmplist.append((float)(part))
				numberlist.append(tmplist)
		returnlist = numberlist
	else:
		returnlist = wordlist	

	return returnlist


#Opens file specified by parameters and returns content as dict or list of dicts in the case of a ROI file
def readfile(source, roi=False, duration=2):
	"""Opens file specified by parameters and returns content as dict or list of dicts in the case of a ROI file"""
	filecontentlist = []
	filecontent = {}
	with open(filename(source, roi, duration), "r") as measurementfile:
		filecontentlist = measurementfile.readlines()
		if (not roi):		
			filecontent.update({"Realtime":(float)(filecontentlist[0].replace(',','').split("\r")[0])})
			filecontent.update({"Livetime":(float)(filecontentlist[1].replace(',','').split("\r")[0])})
			
			for linenumber in xrange(len (filecontentlist)-2):
				filecontent.update({linenumber:(float)(filecontentlist[linenumber + 2].split("\r")[0])})
		else:
			#Set columnnames according to file in readable form given by dict columnnames_roi
			columnnames = []
			for columnname_file in parse_line(filecontentlist[1]):
				columnnames.append(columnnames_roi[columnname_file])
			number_columns = len (columnnames)
			roi = []
			roi_tmp = {}
			for linenumber in xrange(len(filecontentlist) -2):
				roi_tmp = {}
				line_tmp = parse_line(filecontentlist[linenumber+2], True)
				for columnnumber in xrange(len(columnnames)):
					roi_tmp.update({columnnames[columnnumber]:line_tmp[columnnumber]})
				roi.append(roi_tmp)
		
			return roi
	return filecontent


def saveCurrentActivities(date = datetime.date(2017,07,04)):
	strdate = str(date.year) + "-" + str(date.month) + "-" + str(date.day)
	with open("activities_" + strdate + ".txt", "w") as savefile:
		savefile.write(strdate + ";Co;" + str(calc.activity_Co(date)) + "\n")
		savefile.write(strdate + ";Ba;" + str(calc.activity_Ba(date)) + "\n")
		savefile.write(strdate + ";Bi;" + str(calc.activity_Bi(date)) + "\n")

def readCurrentActivities(date = datetime.date(2017,07,04)):
	activities = {}
	strdate = str(date.year) + "-" + str(date.month) + "-" + str(date.day)
	with open("activities_" + strdate + ".txt", "r") as savefile:
		for line in savefile.readlines():
			tmp_line = line.strip("\n").split(";")
			activities.update({tmp_line[1]:float(tmp_line[2])})
	return activities

def loadOutputFile():
	"""data = {"Rn2":[counts, e_counts, [literature_energy_peak, branching_ratio]]}"""
	datei=open("output.txt","r")
	lines=datei.readlines()
	current=""
	data={}
	for i in range(len(lines)):
		    if lines[i].find("Source")!=-1:
		        current = lines[i].split(": ")[1].split("\n")[0]
		        data.update({current:[]})
		    else:
		        points=lines[i].split("\n")[0].split(";")
		        points2=[]
		        for j in points:
		            if j.find("[")==-1 and j.find("]")==-1:
		                points2.append(float(j))
		            elif j.find("[")!=-1:
		                points2.append([float(k) for k in j.split("[")[1].split("]")[0].split(",")])
		        data[current].append(points2)
	datei.close()
	return data

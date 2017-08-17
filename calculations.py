import matplotlib.pyplot as plt
import numpy as np
import math as m
import time
import datetime
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import file_operations as fops
import constants as c
import fit_functions as funcs
import copy

filename_energy_calibration = "parameter_energycalibration.txt"

#Converts the dict given by readfile to plotable list (omitting real and livetime)
def dict_tolist(dict_in):
	"""Converts the dict given by readfile to plotable list (omitting real and livetime)"""
	returnlist = []
	for i in xrange(len(dict_in)-4):
		try:
			returnlist.append(dict_in[i+2])
		except KeyError:
			print "Error at index: " + (str)(i)
	return returnlist

def get_times(source, duration = 2):
	times = {}
	data = fops.readfile(source, False, duration)
	times.update({"Realtime":data["Realtime"]})
	times.update({"Livetime":data["Livetime"]})
	return times

#Returns a list containing 2 lists, with the channel numbers and the energies. Valuepairs are sorted by channelnumber
def energy_calibration(sources = ["Co", "Ba", "Bi"], errors = False):
	"""Returns a list containing 2 lists, with the channel numbers and the energies. Valuepairs are sorted by channelnumber"""
	allroifiles = []
	channel_energy = []
	for source in sources:
		allroifiles.append(fops.readfile(source, True))
		
	for no_roifile in xrange(len(allroifiles)):
		roifile = allroifiles[no_roifile]
		no_peaks = len(roifile)
		fitting_gammas = [[]]
		fitting_gammas = c.gammas[sources[no_roifile]]
		no_gammas = len(fitting_gammas)
		branching_cutoff = 5
		no_popped_gammas = 0
		
		while (no_peaks != no_gammas):
			
			for no_gamma in xrange(no_gammas):
				
				if (fitting_gammas[no_gamma - no_popped_gammas][1] < branching_cutoff):
					fitting_gammas.pop(no_gamma - no_popped_gammas)
					no_popped_gammas += 1
				
			branching_cutoff += 1
			if (branching_cutoff > 100 or no_gammas < no_peaks):
				no_peaks = 0
				no_gammas = 0
				raise Exception("Could not match peaks to gammas.")
				
			no_gammas = len(fitting_gammas)
		
		for no_roi in xrange(len(roifile)):
			tmp_channel_energy = []
			tmp_channel_energy.append(roifile[no_roi]["Centroid"])
			tmp_channel_energy.append(fitting_gammas[no_roi][0])
			if (errors):
				tmp_channel_energy.append(roifile[no_roi]["Uncertainty"])
			channel_energy.append(tmp_channel_energy)
			
		channel_energy.sort()
			
	return channel_energy


def save_energy_calibration(sources= ["Ba", "Co", "Bi"]):
	channel_energy = np.array(energy_calibration(sources,True))
	weighting = []
	for error in channel_energy[:,2]:
		weighting.append(1/error)
	linfit_coeff = np.polyfit(channel_energy[:,0],channel_energy[:,1],1,w=weighting,full=False, cov=True)
	errors = []
	errors.append(m.sqrt(linfit_coeff[1][0][0]))
	errors.append(m.sqrt(linfit_coeff[1][1][1]))
	with open(filename_energy_calibration,"w") as energyfile:
		energyfile.write("Slope;Error Slope;y-Intercept;Error y-Intercept\n")
		energyfile.write((str)(linfit_coeff[0][0]))
		energyfile.write(";")
		energyfile.write((str)(errors[0]))
		energyfile.write(";")
		energyfile.write((str)(linfit_coeff[0][1]))
		energyfile.write(";")
		energyfile.write((str)(errors[1]))
		energyfile.write("\n")
		sourcestmp = "Used Sources;"
		for source in sources:
			sourcestmp += source + ", "
		energyfile.write(sourcestmp.strip(", ") + "\n")


def load_energy_calibration():
	"""Returns the fit parameters of the energycalibration, which were stored in a local file. [slope, err_slope, y-intercept, err_y-intercept]"""
	parameters = []
	with open(filename_energy_calibration,"r") as energyfile:
		filecontent = energyfile.readlines()
		for parameter in filecontent[1].strip("\n").split(";"):
			parameters.append((float)(parameter))
	return parameters

def energy_calibration_fit(sources):
	"""Returns the current fit parameters with errors. [Slope, Error Slope, Y-Intercept, Error Y-Intercept]"""
	save_energy_calibration(sources)
	return load_energy_calibration()

def energy(channel):
	"""Returns energy corresponding to given channel number"""
	fit_param = load_energy_calibration()
	energy_value = fit_param[0]*channel + fit_param[1]
	return energy_value

def channel(energy):
	"""Returns energy corresponding to given channel number"""
	fit_param = load_energy_calibration()
	channel_value = energy / fit_param[0] - fit_param[1]
	return int(channel_value)

def spec_energy(source, duration = 2, substract_background = False):
	"""Returns spectrum of specified source depending on the energy"""
	rawdata = fops.readfile(source, False, duration)
	data = dict_tolist(rawdata)
	spec = []
	backgroundlist = background(rawdata["Livetime"])
	for no_bin in xrange(len(data)):
		tmp_energy = energy(no_bin)
		tmp_data = data[no_bin]
		if (substract_background):
			tmp_data -= backgroundlist[no_bin][1]
		spec.append([tmp_energy, tmp_data])
	return spec

def background_perbin(bin_no, livetime):
	rawfilecontent = fops.readfile("Background")
	return (rawfilecontent[(str)(bin_no)]/rawfilecontent["Livetime"])*livetime

def background(livetime, useenergy = False):
	rawfilecontent = fops.readfile("Background")
	filecontent = dict_tolist(rawfilecontent)
	back_list = []
	for bin_no in xrange(len(filecontent)):
		tmp_x = bin_no
		tmp_y = (filecontent[bin_no]/rawfilecontent["Livetime"])*livetime
		if (useenergy):
			tmp_x = energy(bin_no)
		back_list.append([tmp_x, tmp_y])
	return back_list

def get_double_gaus_starts(x,y,peaks):
    """ calculates starting points for double gaussian"""
    x0 =channel(peaks[0])
    x1=channel(peaks[1])
    mu0 =peaks[2]
    mu1=peaks[3]
    a0 = y[channel(mu0)]
    a1 = y[channel(mu1)]
    sigma0=abs(mu0-mu1)/2.0
    sigma1=abs(mu0-mu1)/2.0
    m=(y[x1]-y[x0])/float(x[x1]-x[x0])
    y0=0
    return [x[x0:x1],y[x0:x1],[a0,mu0,sigma0,a1,mu1,sigma1,m,y0]]
    
def Get_pik_counts(values,e_values):
    
    #geting counts
    counts = abs(values[0] * np.sqrt( 2 * np.pi) *  values[2])
    
    #calculating error
    e_counts =  ( np.sqrt( 2 * np.pi)**2 * ( (e_values[0] * values[2])**2 + (e_values[2] * values[0] )**2))
    e_counts=np.sqrt(e_counts)
    
    #returning Array with 2 indizic.info_Baes 1. number of counts 2. statistic error ( [ number of counts , statistic error ] )
    return [ counts /float(energy(2)-energy(1)), e_counts/float(energy(2)-energy(1)) ]
    
    
def Get_piks(x0,x1,source):
    piks=[]
    for peak in c.gammas[source]:
        if peak[0]-0.7>x0 and peak[0]<x1+0.7:
            piks.append(peak)
    return piks
    

#Peak Counts berechnen durch zaehlen der Counts pro Bin unter der Gauss Kurve von -3sigma bis 3 sigma
def Get_pik_counts_v2(values, e_values, spec_data_list):
	a = values[0]
	mu = values[1]
	sigma = values[2]
	background_gauss = values[3]
	
	bin_from = channel(mu - 3*sigma)
	bin_to = channel(mu + 3*sigma)
	
	counts = 0
	
	for bin_no in xrange(bin_from, bin_to):
		counts += spec_data_list[bin_no] - background_gauss
	
	print "Peak Counts v2 at Peak " + (str)(mu) + ": " + (str)(counts)
	return [counts,e_values]

    
def calc_energyefficiency_points():
    datei = open("output.txt","r")
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
                        points2.append([float(k) for k in j.split("[")[1].split("]")[0].split(", ")])
                data[current].append(points2)
    x=[]
    y=[]
    x2=[[],[],[]]
    y2=[[],[],[]]
    error=[]
    for source in ["Ba","Co","Bi"]:
        for peak in data[source]:
            counts=peak[0]
            errorcounts=peak[1]
            branching=0
            energie=0
            if source=="Ba":
                Activity=c.info_Ba["Activity"]*1000
                lifetime=c.info_Ba["Lifetime"]
                timestr=c.info_Ba["ActTime"].split(".")
                time1 = time.mktime(datetime.datetime.strptime(timestr[0]+"/"+timestr[1]+"/"+timestr[2], "%d/%m/%Y").timetuple())
                time2 = time.mktime(datetime.datetime.strptime("04/07/2017", "%d/%m/%Y").timetuple())
                dt=time2-time1
                halftimes= dt/float(c.info_Ba["Halflife"]*3600*24*365.25)
                percentage = 0.5**halftimes
                angel =m.sin(m.atan(3.5)/2)**2
            if source=="Co":
                Activity=c.info_Co["Activity"]*1000
                timestr=c.info_Co["ActTime"].split(".")
                lifetime=c.info_Co["Lifetime"]
                time1 = time.mktime(datetime.datetime.strptime(timestr[0]+"/"+timestr[1]+"/"+timestr[2], "%d/%m/%Y").timetuple())
                time2 = time.mktime(datetime.datetime.strptime("04/07/2017", "%d/%m/%Y").timetuple())
                dt=time2-time1
                halftimes= dt/float(c.info_Co["Halflife"]*3600*24*365.25)
                percentage = 0.5**halftimes
                angel =m.sin(m.atan(3.5)/2)**2
            if source=="Bi":
                Activity=c.info_Bi["Activity"]*1000
                timestr=c.info_Bi["ActTime"].split(".")
                lifetime=c.info_Bi["Lifetime"]
                time1 = time.mktime(datetime.datetime.strptime(timestr[0]+"/"+timestr[1]+"/"+timestr[2], "%d/%m/%Y").timetuple())
                time2 = time.mktime(datetime.datetime.strptime("04/07/2017", "%d/%m/%Y").timetuple())
                dt=time2-time1
                halftimes= dt/float(c.info_Bi["Halflife"]*3600*24*365.25)
                percentage = 0.5**halftimes
                angel =m.sin(m.atan(3.5)/2)**2
                
            for singlepeak in range(len(peak)-2):
                branching+= peak[singlepeak+2][1]/100.0
            for singlepeak in range(len(peak)-2):
                energie+= peak[singlepeak+2][1]*peak[singlepeak+2][0]/100.0/branching
            Activity=Activity*percentage*angel*branching*lifetime
            print energie,Activity,counts,branching,percentage,m.atan(3.5)*180/m.pi
            x.append(energie)
            y.append(counts/Activity)
            error.append(m.sqrt((errorcounts/Activity)**2+(counts/Activity*0.15)**2))
            if source=="Bi":
                x2[0].append(energie)
                y2[0].append(counts/Activity)
            if source=="Ba":
                x2[1].append(energie)
                y2[1].append(counts/Activity)
            if source=="Co":
                x2[2].append(energie)
                y2[2].append(counts/Activity)
            
    a=True
    while a:
        a=False
        for i in range(len(x)-1):
                if x[i+1]<x[i]:
                    a=True
                    a=x[i]
                    b=y[i]
                    cc=error[i]
                    x[i]=x[i+1]
                    y[i]=y[i+1]
                    error[i]=error[i+1]
                    x[i+1]=a
                    y[i+1]=b
                    error[i+1]=cc
    for i in range(len(x)):
        print x[i],y[i],error[i]
    return [x,y,error,x2,y2]


def efficienzy(energie):
	#x,y ,err ,a,b= calc_energyefficiency_points()
	data = np.array(efficiency_perpeak())
	data = data[data[:, 0].argsort()]
	x = data[:,0]
	y = data[:,1]
	y_err = data[:,2]
 
	n = 0
	
	if (energie > x[-1]):
		return 0.0,0.0
	
	
	while x[n] < energie:
		#print energie
		#print x[n]
		n += 1
	m = (y[n]-y[n-1]) / (x[n]-x[n-1])
	t = y[n] - m*x[n]
	m_err= (y_err[n] / (x[n]-x[n-1]))**2
	m_err+=  (y_err[n-1] / (x[n]-x[n-1]))**2
	m_err=np.sqrt(m_err)
	t_err = y_err[n]**2
	t_err += (m_err*x[n])**2
	t_err =np.sqrt(t_err)
	total_err=(m_err*energie)**2+t_err**2   
 
     
	if (energie <= x[0]):
		return 0.0,0.0
	return m*energie + t, np.sqrt(total_err)
 
    
    
def get_actual_decays(counts,e_counts,energie,branching):
    data=open( "Raumwinkel.txt","r")
    a= data.readlines()
    winkel=float(a[0].split("\n")[0])
    e_winkel=float(a[1])
    eff , e_eff = efficienzy(energie)
    
    count=counts/eff*m.pi*4/winkel/branching
    e_count_angle=counts/branching/eff*m.pi*4/winkel/winkel*e_winkel
    e_count_eff=counts/branching/eff*m.pi*4/winkel/eff*e_eff
    e_count_counts=e_counts/branching/eff*m.pi*4/winkel
    
    return [count,m.sqrt((e_count_angle)**2+(e_count_eff)**2+e_count_counts**2)]
   

def activity_Co(date):
	"""Calculates Activity in kBq at given date."""
	ref_date = datetime.datetime.strptime(c.info_Co["ActTime"],"%d.%m.%Y").date()
	halflife = c.info_Co["Halflife"]*365
	ref_activity = c.info_Co["Activity"]
	diff = abs((date - ref_date).days)
	decay_constant = m.log(2)/halflife
	
	activity = ref_activity * m.exp(-decay_constant*diff)
	
	return activity

def activity_Ba(date):
	"""Calculates Activity in kBq at given date."""
	ref_date = datetime.datetime.strptime(c.info_Ba["ActTime"],"%d.%m.%Y").date()
	halflife = c.info_Ba["Halflife"]*365
	ref_activity = c.info_Ba["Activity"]
	diff = abs((date - ref_date).days)
	decay_constant = m.log(2)/halflife
	
	activity = ref_activity * m.exp(-decay_constant*diff)
	
	return activity

def activity_Bi(date):
	"""Calculates Activity in kBq at given date."""
	ref_date = datetime.datetime.strptime(c.info_Bi["ActTime"],"%d.%m.%Y").date()
	halflife = c.info_Bi["Halflife"]*365
	ref_activity = c.info_Bi["Activity"]
	diff = abs((date - ref_date).days)
	decay_constant = m.log(2)/halflife
	
	activity = ref_activity * m.exp(-decay_constant*diff)
	
	return activity


def efficiency_perpeak():
    
    
	"""Returns literature values of [energy of peak, efficiency of peak]"""
	raumwinkel_calibration = 0.36
	sources = ["Co", "Bi", "Ba"]
	activities = fops.readCurrentActivities()
	peaks = []
	efficiency = []
	for source in sources:
		livetime = get_times(source)["Livetime"]
		fitted_peaks = fops.loadOutputFile()[source]				#[[counts, e_counts, [energy, branching]], counts2, e_counts2, [...,...] ,...]
		tmp_activity = activities[source]*1000							#kBq -> Bq = 1/s
		tmp_gammas = c.gammas[source]
		
		for no_gamma in xrange(len(fitted_peaks)):
			gamma = fitted_peaks[no_gamma]
			
			tmp_branching = gamma[2][1]/100.0
			for i in range(len(gamma)-3):							#Double Peak?
				tmp_branching += gamma[i+3][1]/100.0
			
			
			tmp_act_pergamma = tmp_activity * tmp_branching			#activity * branching ratio(in percent)/100
			peaks.append([gamma[2][0], tmp_act_pergamma])			#[energy of gamma, activity per gamma]
			
			tmp_counts_persecond = gamma[0]/livetime				#counts per second = counts measured / livetime
			tmp_efficiency = tmp_counts_persecond / (tmp_act_pergamma * raumwinkel_calibration)			#counts per second / activity per gamma
			tmp_efficiency_err=gamma[1]/livetime/ (tmp_act_pergamma * raumwinkel_calibration)
			efficiency.append([gamma[2][0], tmp_efficiency, tmp_efficiency_err])		#[energy of gamma, efficiency]
			
			#print "energy:" + str(gamma[2][0])
			#print "act: " + str(tmp_activity)
			#print "branch: " + str(gamma[2][1])
			#print "act_pergamma: " + str(tmp_act_pergamma)
			#print "counts_per second: " + str(tmp_counts_persecond)
			#print "efficiency: " + str(tmp_efficiency)
			
			
	return efficiency;




def Pb_calculator():
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
	outlines=["Source : Rn2\n"]
	for peak in data["Rn2"]:
		counts , e_counts=get_actual_decays(peak[0],peak[1],peak[2][0],peak[2][1]/100.0)
		outlines.append(str(counts)+";"+str(e_counts)+";["+str(peak[2][0])+","+str(peak[2][1])+"]\n")
	outlines.append("Source : Rn4\n")
	for peak in data["Rn4"]:
		counts , e_counts=get_actual_decays(peak[0],peak[1],peak[2][0],peak[2][1]/100.0)
		outlines.append(str(counts)+";"+str(e_counts)+";["+str(peak[2][0])+","+str(peak[2][1])+"]\n")
	datei=open("totaldecays.txt","w")
	datei.writelines(outlines)
	datei.close()




def get_Pb_conzentration(counts,epsf,V,t1,t2,T,lambdaPo,lambdaPb,e_counts):
	#aufbroeseln der Gleichung in eigene Klammern
    	print T,lambdaPo,lambdaPb
     
     
     
	part1=counts/epsf/V
	part2=1/lambdaPo+1/lambdaPb+lambdaPb/(lambdaPo*(lambdaPo-lambdaPb))+lambdaPo*m.exp(-lambdaPb*T)/(lambdaPb*(lambdaPb-lambdaPo))
	part3=m.exp(-lambdaPb*t1)-m.exp(-lambdaPb*t2)
	part4=lambdaPb/(lambdaPo*(lambdaPb-lambdaPo))*(1-m.exp(-lambdaPo*T))
	part5=m.exp(-lambdaPo*t1)-m.exp(-lambdaPo*t2)

	error=part1/(part2*part3+part4*part5)*e_counts/counts
	return [part1/(part2*part3+part4*part5),error]



def get_Pb14_konz(counts,e_counts):

	lambdaPo=c.info_Po_218["lambda"]
	lambdaPb=c.info_Pb_214["lambda"]
	t1=c.info_Rn2["Starttime"]
	t2=c.info_Rn2["Endtime"]
	V=c.info_Filter["V"]
	epsif=c.info_Filter["epsi"]
	T=c.info_Filter["Filtertime"]
	konzentration=[]
	e_konzentration=[]
	for count in range(len(counts)):
		konz,e_konz=get_Pb_conzentration(counts[count],epsif,V,t1,t2,T,lambdaPo,lambdaPb,e_counts[count])
		konzentration.append(konz)
		e_konzentration.append(e_konz)
	return [konzentration , e_konzentration]


def get_Pb12_konz(counts,e_counts):

	lambdaPo=c.info_Po_216["lambda"]
	lambdaPb=c.info_Pb_212["lambda"]
	t1=c.info_Rn4["Starttime"]
	t2=c.info_Rn4["Endtime"]
	V=c.info_Filter["V"]
	epsif=c.info_Filter["epsi"]
	T=c.info_Filter["Filtertime"]
	konzentration=[]
	e_konzentration=[]
	for count in range(len(counts)):
		konz,e_konz=get_Pb_conzentration(counts[count],epsif,V,t1,t2,T, lambdaPo,lambdaPb,e_counts[count])
		konzentration.append(konz)
		e_konzentration.append(e_konz)
	return [konzentration , e_konzentration]

def save_Pb_konzentrations():
	datei=open("totaldecays.txt","r")
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
	
	countsPb14=[]
	e_countsPb14=[]
	countsPb12=[]
	e_countsPb12=[]
	for peak in data["Rn2"]:
		countsPb14.append(peak[0])
		e_countsPb14.append(peak[1])
	for peak in data["Rn4"]:
		countsPb12.append(peak[0])
		e_countsPb12.append(peak[1])


	Pb14 , e_Pb14=get_Pb14_konz(countsPb14,e_countsPb14)
	Pb12 , e_Pb12=get_Pb12_konz(countsPb12,e_countsPb12)
	outlines=["Source : Pb14\n"]

	for i in range(len(Pb14)):
		outlines.append(str(Pb14[i])+";"+str(e_Pb14[i])+"\n")
	outlines.append("Source : Pb12\n")

	for i in range(len(Pb12)):
		outlines.append(str(Pb12[i])+";"+str(e_Pb12[i])+"\n")

	datei=open("Pbkonzentration.txt","w")
	datei.writelines(outlines)
	datei.close()
	
def get_Radon_konzentration(lambdaRn,lambdaPb,conzPb,e_conzPb):
	return [conzPb*lambdaPb/lambdaRn,e_conzPb*lambdaPb/lambdaRn]


def save_Radon_conzentration():
	datei=open("Pbkonzentration.txt","r")
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

	outlines=["Source : Rn222\n"]
	for peak in data["Pb14"]:
		count,e_count=get_Radon_konzentration(c.info_Rn_222["lambda"],c.info_Pb_214["lambda"],peak[0],peak[1])
		outlines.append(str(count)+";"+str(e_count)+"\n")
	outlines.append("Source : Rn220\n")
	for peak in data["Pb12"]:
		count,e_count=get_Radon_konzentration(c.info_Rn_220["lambda"],c.info_Pb_212["lambda"],peak[0],peak[1])
		outlines.append(str(count)+";"+str(e_count)+"\n")
	datei=open("Rnkonzentration.txt","w")
	datei.writelines(outlines)
	datei.close()
 
 
def get_chiquad(x,y,funktion,parameter,sigma):
    chi=0
    for i in range(len(x)):
        if len(parameter)==0:
            chi+=(funktion(x[i])-y[i])**2
        elif len(parameter)==1:
            chi+=(funktion(x[i],parameter[0])-y[i])**2
        elif len(parameter)==2:
            chi+=(funktion(x[i],parameter[0],parameter[1])-y[i])**2
        elif len(parameter)==3:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2])-y[i])**2
        elif len(parameter)==4:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2],parameter[3])-y[i])**2
        elif len(parameter)==5:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2],parameter[3],parameter[4])-y[i])**2
        elif len(parameter)==6:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2],parameter[3],parameter[4],parameter[5])-y[i])**2
        elif len(parameter)==7:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2],parameter[3],parameter[4],parameter[5],parameter[6])-y[i])**2
        elif len(parameter)==8:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2],parameter[3],parameter[4],parameter[5],parameter[6],parameter[7])-y[i])**2
        elif len(parameter)==9:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2],parameter[3],parameter[4],parameter[5],parameter[6],parameter[7],parameter[8])-y[i])**2
        elif len(parameter)==10:
            chi+=(funktion(x[i],parameter[0],parameter[1],parameter[2],parameter[3],parameter[4],parameter[5],parameter[6],parameter[7],parameter[8],parameter[9])-y[i])**2
    return chi/(sigma*sigma)
    
def get_parameter_error(x,y,funktion,parameter,sigma,steps=[0.0001]):
    errors=[]
    st=steps
    if len(steps)==1:
        st=[]
        for i in range(len(parameter)):
            st.append(steps[0])
    chi0=get_chiquad(x,y,funktion,parameter,sigma)
    print "chiquad=",chi0,len(x)
    for i in range(len(parameter)):
        
        errors.append([])
        para=copy.deepcopy(parameter)
        while abs(get_chiquad(x,y,funktion,para,sigma)-chi0)<1.0:
            para[i]-=st[i]
        errors[-1].append(copy.copy(para[i]))
        errors[-1].append(copy.copy(parameter[i]))
        para=copy.deepcopy(parameter)
        while abs(get_chiquad(x,y,funktion,para,sigma)-chi0)<1.0:
            para[i]+=st[i]
        errors[-1].append(copy.copy(para[i]))
    return errors
        
    
    
import matplotlib.pyplot as plt
import matplotlib.lines as legendlines
import numpy as np

import datetime

from matplotlib.backends.backend_pdf import PdfPages

import file_operations as fops
import constants as c
import calculations as calc
import fit_functions as funcs
import raum

printable_sources = {"Ba":"$^{133}$Ba", "Co":"$^{60}$Co", "Bi":"$^{207}$Bi", "Rn":"Filterpaper", "Background":"Background"}



cutoff_energy = 2000				#keV


def plot_energycalibration(redo_fit=False, sources=["Co", "Ba", "Bi"]):
	
	#redo the fit?
	if (redo_fit):
		calc.save_energy_calibration(sources)
	
	fit_params = calc.load_energy_calibration()
	
	calib_points = np.array(calc.energy_calibration(sources))
	
	plt.scatter(calib_points[:,0], calib_points[:,1], marker="x", color="r")
	plt.plot(calib_points[:,0], fit_params[0] * calib_points[:,0] + fit_params[2], "-b")
	
	sources_string = ""
	for source in sources:
		sources_string += "${}^{" + str(c.info[source]["Massnumber"]) + "}$" + source + ", "
	sources_string = sources_string.strip(", ")
	
	#Labels
	title = "Energycalibration\n using " + sources_string
	x_axis = "ADC Channel [1]"
	y_axis = "Energy [keV]"
	
	plt.title(title)
	
	plt.xlabel(x_axis)
	plt.ylabel(y_axis)
	
	plt.show()


def plot_energy_specs_meas_back_noback():
	pages = PdfPages("specs_meas-back-noback.pdf")
	title = ""
	source = ""
	duration = 0
	#for sourceduration in ["Rn2", "Rn4", "Co", "Bi", "Ba"]:
	for sourceduration in ["Rn4"]:
		if (sourceduration.startswith("Rn")):
			source = sourceduration[:2]
			duration = int(sourceduration[-1])
		else:
			source = sourceduration
			
		data = np.array(calc.spec_energy(source, duration, False))
		data_noback = np.array(calc.spec_energy(source, duration, True))
		data_back = np.array(calc.background(calc.get_times(source)["Livetime"], True))

		x = data[:calc.channel(cutoff_energy),0]
		y = data[:calc.channel(cutoff_energy),1]

		x_nb = data_noback[:calc.channel(cutoff_energy),0]
		y_nb = data_noback[:calc.channel(cutoff_energy),1]

		x_b = data_back[:calc.channel(cutoff_energy),0]
		y_b = data_back[:calc.channel(cutoff_energy),1]

		plt.plot(x, y, "-b")
		plt.plot(x_b, y_b, ":g")
		plt.plot(x_nb, y_nb, "-r")
		plt.show()
		if (source == "Rn"):
			title = printable_sources[source] + ", " + (str)(duration) + "h, Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
		else:
			title = printable_sources[source] + ", Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
		plt.title(title)
		plt.xlabel("Energy [keV]")
		plt.ylabel("Counts")
		
		measurement_line = legendlines.Line2D([],[], color="blue", label=printable_sources[source])
		background_line = legendlines.Line2D([],[], color="green", label="Background")
		measurement_noback_line = legendlines.Line2D([],[], color="red", label=printable_sources[source] + " w/o Bkg")
		
		plt.legend(handles=[measurement_line, background_line, measurement_noback_line])
		
		pages.savefig()
		#funcs.Test_gaus_plot(x_nb[calc.channel(230):calc.channel(241)],y_nb[calc.channel(230):calc.channel(241)])
		#funcs.Test_double_gaus_plot(x_nb[calc.channel(61):calc.channel(79)],y_nb[calc.channel(61):calc.channel(79)],71,73,y_nb[calc.channel(71)]-y_nb[calc.channel(73)],y_nb[calc.channel(61)]-y_nb[calc.channel(79)])
	pages.close()
	
	
def plot_energy_specs_noback(plot_hist=False, channels_per_bin=10):
	pages = PdfPages("specs_noback.pdf")
	title = ""
	source = ""
	duration = 0
	for sourceduration in ["Rn2", "Rn4", "Co", "Bi", "Ba"]:
		if (sourceduration.startswith("Rn")):
			source = sourceduration[:2]
			duration = int(sourceduration[-1])
		else:
			source = sourceduration
		
		if (plot_hist):			#Plot as Histogram
			
			data_noback = np.array(calc.spec_energy(source, duration, True))[:calc.channel(cutoff_energy),:]
			
			hist_bins = len(data_noback)
			
			data_hist = hist_values(data_noback)
			
			plt.hist(data_hist, hist_bins/channels_per_bin, facecolor="red")
			
		else:					#Plot as Function
			
			data_noback = np.array(calc.spec_energy(source, duration, True))
			
			x_nb = data_noback[:calc.channel(cutoff_energy),0]
			y_nb = data_noback[:calc.channel(cutoff_energy),1]

			plt.plot(x_nb, y_nb, "-r")
		
		if (source == "Rn"):
			title = printable_sources[source] + ", " + (str)(duration) + "h, Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
		else:
			title = printable_sources[source] + ", Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
		plt.title(title)
		plt.xlabel("Energy [keV]")
		plt.ylabel("Counts")
		
		measurement_noback_line = legendlines.Line2D([],[], color="red", label=printable_sources[source] + " w/o Bkg")
		
		plt.legend(handles=[measurement_noback_line])
		
		pages.savefig()
		plt.close()
          
	pages.close()


def peakplotter(show_plots = True, plot_hist=False):
	
	pages = PdfPages("calibration_peak_plots.pdf")
	title = ""
	source = ""
	duration = 0
	output = open("output.txt","w")
     	outputlines=[]
	for sourceduration in ["Rn2","Rn4","Co", "Bi", "Ba"]:   
         
		if (sourceduration.startswith("Rn")):
			source = sourceduration[:2]
			duration = int(sourceduration[-1])
		else:
			source = sourceduration
		#saving sourcename in outputfile
    		outputlines.append("Source : "+sourceduration+"\n")     
			
		data_noback = np.array(calc.spec_energy(source, duration, True))
		
		x_nb = data_noback[:calc.channel(cutoff_energy),0]
		y_nb = data_noback[:calc.channel(cutoff_energy),1]


            #Plotting peaks with only one peak in range
		for peak in c.calibration_peaks_single_range[sourceduration]:
                	x_pk = x_nb[calc.channel(peak[0]):calc.channel(peak[1])]
                	y_pk = y_nb[calc.channel(peak[0]):calc.channel(peak[1])]
                	values , errors = funcs.Fit_Gaussian(x_pk,y_pk,background=True)
                	sigma=0
                	for i in range(len(x_pk)):
                         sigma+=(funcs.Gauss(x_pk[i], values[0],values[1],values[2],values[3])-y_pk[i])**2
                	sigma=np.sqrt(sigma/float(len(x_pk)-1))
                	error2= calc.get_parameter_error(x_pk,y_pk,funcs.Gauss,values,sigma)
                	err=[]
                	err2=[]
                	for i in range(len(error2)):
                         err.append((error2[i][1]-error2[i][0])/2.0)
                	print errors
                	print err
                	print err2
                	print error2
                	print "\n\n\n\n\n"
                 
                 
                	x_gaus=[]
                	y_gaus=[]
                	for i in range(1001):
                            x_gaus.append(x_pk[0]+i*(x_pk[-1]-x_pk[0])/1000.0)
                            y_gaus.append( funcs.Gauss(x_gaus[-1], values[0],values[1],values[2],values[3]))
    
    
                
                 	ax=plt.figure("test")
                 	
                 	#Option to plot spectrum as histogram
                 	if (plot_hist):
						
						hist_bins = len(x_pk)
						
						xypairs = map(list, zip(x_pk, y_pk)	)			#[x0,x1,x2,...];[y0,y1,y2,...] to [[x0,y0],[x1,y1],[x2,y2],...]
						
						hist_val = hist_values(xypairs)
						
						plt.hist(hist_val, hist_bins, facecolor="blue")
                 	
                 	else:
						plt.plot(x_pk, y_pk, "-b")
                	
                	plt.plot(x_gaus, y_gaus, "-r")
    		
                	if (source == "Rn"):
            			title = printable_sources[source] + ", " + (str)(duration) + "h, Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
                	else:
            			title = printable_sources[source] + ", Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
                	plt.title(title)
                	plt.xlabel("Energy [keV]")
                	plt.ylabel("Counts")
                	textstr = '$\mathrm{a =}%.2f \pm %.2f$\n$\mu=%.2f \pm %.2f$\n$\sigma=%.2f \pm %.2f$\n$\mathrm{y0 =}%.2f \pm %.2f$'%(values[0],errors[0], values[1],errors[1], values[2],errors[2], values[3],errors[3])
                	ax.text(0.17, 0.85, textstr, fontsize=14,verticalalignment='top',bbox={"facecolor":"white","alpha":0.5,"pad":10})
                 
                	measurement_noback_line = legendlines.Line2D([],[], color="red", label=printable_sources[source] + " w/o Bkg")
                 
                	plt.legend(handles=[measurement_noback_line])
                	pages.savefig()
                	if (show_plots):
						plt.show()
                	plt.close()
                 
                 #saving data in output file
                	piks = calc.Get_piks(x_pk[0],x_pk[-1],sourceduration)
                	#print piks,x_pk[0],x_pk[-1]
                	counts ,e_counts = calc.Get_pik_counts(values[0:3],err[0:3])
                	outputlines.append(str(counts) +";"+str(e_counts))
                	for peak in piks:
						outputlines[-1]+= ";"+str(peak)
                	outputlines[-1]+="\n"                    
                 
                 
              
            #Plotting peaks with 2 peaks in range
		for peak in c.calibration_peaks_double_range[sourceduration]:
      
      
                	x_pk , y_pk , start = calc.get_double_gaus_starts(x_nb,y_nb,peak)
                 
                	values , errors = funcs.Fit_double_Gaussian(x_pk,y_pk,start)
                 
                	sigma=0
                	for i in range(len(x_pk)):
                         sigma+=(funcs.double_Gauss(x_pk[i], values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7])-y_pk[i])**2
                	sigma=np.sqrt(sigma/float(len(x_pk)-1))
                	error2= calc.get_parameter_error(x_pk,y_pk,funcs.double_Gauss,values,sigma)
                	err=[]
                	err2=[]
                	for i in range(len(error2)):
                         err.append((error2[i][1]-error2[i][0])/2.0)
                	print errors
                	print err
                	print err2
                	print error2
                	print "\n\n\n\n\n"
                 
                 
                	x_gaus=[]
                	y_gaus=[]
                	for i in range(1001):
                            x_gaus.append(x_pk[0]+i*(x_pk[-1]-x_pk[0])/1000.0)
                            y_gaus.append( funcs.double_Gauss(x_gaus[-1], values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7]))
    
    
                
                 	ax=plt.figure("test")
                 	
                 	
                 	#Option to plot spectrum as histogram
                 	if (plot_hist):
						
						hist_bins = len(x_pk)
						
						xypairs = map(list, zip(x_pk, y_pk)	)			#[x0,x1,x2,...];[y0,y1,y2,...] to [[x0,y0],[x1,y1],[x2,y2],...]
						
						hist_val = hist_values(xypairs)
						
						plt.hist(hist_val, hist_bins, facecolor="blue")
                 	
                 	else:
						plt.plot(x_pk, y_pk, "-b")
						
                	plt.plot(x_gaus, y_gaus, "-r")
    		
                	if (source == "Rn"):
            			title = printable_sources[source] + ", " + (str)(duration) + "h, Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
                	else:
            			title = printable_sources[source] + ", Livetime " + (str)(calc.get_times(source, duration)["Livetime"]) + ", Realtime " + (str)(calc.get_times(source, duration)["Realtime"])
                	plt.title(title)
                	plt.xlabel("Energy [keV]")
                	plt.ylabel("Counts")
                	textstr = '$\mathrm{a =}%.2f $\n$\mu=%.2f $\n$\sigma=%.2f $\n$\mathrm{a2 =}%.2f $\n$\mu 2=%.2f $\n$\sigma 2=%.2f $\n$\mathrm{m =}%.2f$\n$\mathrm{y0 =}%.2f$'%(values[0], values[1], values[2], values[3],values[4], values[5], values[6], values[7])
                	ax.text(0.17, 0.85, textstr, fontsize=14,verticalalignment='top',bbox={"facecolor":"white","alpha":0.5,"pad":10})
                 
                	measurement_noback_line = legendlines.Line2D([],[], color="red", label=printable_sources[source] + " w/o Bkg")
                 
                	plt.legend(handles=[measurement_noback_line])
                	pages.savefig()
                	#plt.show()
                	plt.close()
                 #saving data in output file
                 
                 
                	mitte= calc.channel((values[1]+values[4])/2.0)
                	piks = calc.Get_piks(x_pk[0],x_pk[np.where(y_pk == y_nb[mitte])[0]],sourceduration)
                 
                 
                	counts ,e_counts = calc.Get_pik_counts(values[0:3],err[0:3])
                	
                	#counts ,e_counts = calc.Get_pik_counts_v2(values,errors, y_nb)
                	
                	if len(piks)!=0:
                        	outputlines.append(str(counts) +";"+str(e_counts))
                        	for peak in piks:
                                 outputlines[-1]+= ";"+str(peak)
                        	outputlines[-1]+="\n"                    
                	piks = calc.Get_piks(x_pk[np.where(y_pk == y_nb[mitte])[0]],x_pk[-1],sourceduration)
                 
                 
                	counts ,e_counts = calc.Get_pik_counts(values[3:6],err[3:6])
                	if len(piks)!=0:
                        	outputlines.append(str(counts) +";"+str(e_counts))
                        	for peak in piks:
                                 outputlines[-1]+= ";"+str(peak)
                        	outputlines[-1]+="\n"                    
	output.writelines(outputlines)
	output.close()            
	pages.close()
	
	
	
def plot_energyefficiency_david():
	data = np.array(calc.calc_energyefficiency_points())
	plt.plot(data[:,0], data[:,1], marker='x')
	plt.show()
	print data

def plot_energyefficiency():
	data = np.array(calc.efficiency_perpeak())
	data = data[data[:, 0].argsort()]
	plt.plot(data[:,0], data[:,1], marker='x')
	#plt.show()
	print data
	with PdfPages("energyefficiency.pdf") as pages:
		pages.savefig()

def plot_efficiency_v3():
	x = []
	y = []
	x0 = 0
	xmax = 2000
	for i in range(1001):
		x.append(x0 + i*(xmax - x0)/1000.0)
		y.append(calc.efficienzy(x[-1])[0])
	
	plt.plot(x, y)
	plt.show()


def hist_values(spec_xy):
	
	"""Input: [[x0, y0], [x1, y1], ...]
	Output: [x0, x0, x0, ..., x1, x1, ...] y0 times x0, y1 times x1, ...
	Converts xy point pairs suitable for pyplot.plot() to data for pyplot.hist()"""
	
	xy = np.array(spec_xy)
	ret = []
	for xypair in xy:
		if (xypair[1] > 0.0):
			for event in xrange(int(xypair[1])):
				ret.append(xypair[0])
	return ret

def plot_hist():
	spec = calc.spec_energy("Co")
	x = hist_values(spec)
	print x
	plt.hist(x,5000)
	plt.show()
plot_energy_specs_meas_back_noback()
raum.do_everything(template=False)
peakplotter()
calc.Pb_calculator()
calc.save_Pb_konzentrations()
calc.save_Radon_conzentration()

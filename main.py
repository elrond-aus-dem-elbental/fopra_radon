#!/usr/bin/python

import calculations as calc
import plot


def run_calc_radon():
	calc.Pb_calculator()
	calc.save_Pb_konzentrations()
	calc.save_Radon_conzentration()


#plot.plot_hist()


#plot.peakplotter(show_plots=False, plot_hist=True)
#plot.plot_energy_specs_noback(plot_hist=True)
plot.plot_energycalibration(True)

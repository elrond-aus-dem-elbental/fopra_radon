#Emitted gammas of the source: [energy in keV, branching_ratio in percent]
#Source: Table of Isotopes - The Lund/LBNL, Nuclear Data Search, Version 2.0, February 1999, S.Y.F. Chu. L.P. Ekstroem and R.B. Firestone
import math as m

gammas = {
	"Ba":[
		[53.161,2.199],
		[79.6139,2.62],
		[80.9971,34.06],
		[276.398,7.164],
		[302.853,18.33],
		[356.017,62.05],
		[383.851,8.94]
		],
	"Co":[
		[1173.237,99.9736],
		[1332.501,99.9856]
		],
	"Bi":[
		[72.805,21.3],
		[74.969,35.8],
		[84.450,4.26],
		[84.938,8.19],
		[87.300,2.98],
		[569.702,97.74],
		[1063.662,74.5],
		[1770.237,6.87]
		],
    "Rn4":[
            [238.6,43.6],
            [300.1,3.34],
            [77.1,17.5]
    ],
    "Rn2":[
            [77.1,10.8],
            [241.9,7.46],
            [295.2,19.2],
            [351.9,37.1],
            [785.9,1.09]
    ]
}



#Plotgrenzen der verschiedenen Peaks
calibration_peaks_single_range= {
	"Ba":[
		[70,90],
		[260,290],
		[290,315],
		[340,370],
		[370,400]
		],
	"Co":[
		[1165,1185],
		[1320,1350]
		],
	"Bi":[
		[1050,1080],
		[550,590],
		[1760,1790]
		],
	"Rn2":[
		[285,310],
		[330,370],
		[740,800]
		],
	"Rn4":[
		[280,310],
		[50,100]
		]
  
}

#Plotgrenzen der Doppelpeaks und Startwerte fuer die Mittelwerte (um Konvergenz sicherzustellen)
calibration_peaks_double_range= {
	"Ba":[
		],
	"Co":[
		],
	"Bi":[
		[61,79,71.5,73.5],
		[79,92,84,87]
		],
	"Rn2":[
		[220,260,238,241]
		],
	"Rn4":[
		[220,260,238,241]
		]
}


info = {"Co": {"Massnumber":60}, "Ba": {"Massnumber":133}, "Bi": {"Massnumber":207}}


#Activities of the source at specified date in kBq, Halflife in years (?) (y 5 <-??? and time of exposer(lifetime))
info_Filter = {"Filtertime":60*85,"V":200/3600.0,"epsi":0.9997}
info_Rn2 = {"Starttime":3*60,"Endtime":6874.158+180.0,"Lifetime":6468.992}
info_Rn4 = {"Starttime":3*60,"Endtime":3*60+14660.,"Lifetime":13786.}
info_Rn_220 = {"Halftime":55.6,"lambda":-m.log(0.5)/55.6}
info_Rn_222 = {"Halftime":3.82*24*3600.,"lambda":-m.log(0.5)/(3.82*24*3600.)}
info_Po_216 = {"Halftime":0.15,"lambda":-m.log(0.5)/0.15}
info_Po_218 = {"Halftime":3.11*60,"lambda":-m.log(0.5)/(3.11*60.)}
info_Pb_212 = {"Halftime":10.64*3600,"lambda":-m.log(0.5)/(10.64*3600.)}
info_Pb_214 = {"Halftime":26.8*60,"lambda":-m.log(0.5)/(26.8*60.)}
info_Ba = {"Activity":3.08, "ActTime":"12.01.2005", "Halflife":10.51,"Lifetime":650.0}
info_Co = {"Activity":47.3, "ActTime":"19.02.1999", "Halflife":5.2714,"Lifetime":151.0}		#Activity was corrected on paper (used to be 4.73) (?)
info_Bi = {"Activity":37, "ActTime":"24.11.2008","Halflife":31.55,"Lifetime":402.6}

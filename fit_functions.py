import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def xponential_decay(x,a,b,c):
    return a* np.exp(-b*x)+c
    
def x_negv_pow(x,a,b):
    return a/(x+b)


def Gauss(x, a, mu, sigma,y0):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2)) + y0
    
    
def double_Gauss(x, a1, mu1, sigma1, a2, mu2, sigma2,m,y0):
    return a1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) + a2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2)) + y0 + m*x
    
def Gauss_nobackground(x, a, mu, sigma):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2))
    
def Fit_exponential_decay(x,y,start):
    x = np.asarray(x)
    y = np.asarray(y)
    #guessing starting point
    #calculating actual fit
    popt,pcov = curve_fit(xponential_decay, x, y, p0=start)
    #recieveing Errors
    #print popt, pcov
    try:
        error=np.sqrt(np.diag(pcov))
    except:
        error=[]
        
    #returns array containing two arrays 1. values 2. corresponding error ( [ [ a , mu , sigma , y0 ] , [ e_a, e_mu , e_ sigma , e_y0 ] ] )
    return [popt,error]
    
def Fit_x_negv_pow(x,y,start):
    x = np.asarray(x)
    y = np.asarray(y)
    #guessing starting point
    #calculating actual fit
    popt,pcov = curve_fit(x_negv_pow, x, y, p0=start)
    #recieveing Errors
    #print popt, pcov
    try:
        error=np.sqrt(np.diag(pcov))
    except:
        error=[]
        
    #returns array containing two arrays 1. values 2. corresponding error ( [ [ a , mu , sigma , y0 ] , [ e_a, e_mu , e_ sigma , e_y0 ] ] )
    return [popt,error]
    
    
    
def Fit_double_Gaussian (x,y,start):
    
    #C reates gaussian Fit
    x = np.asarray(x)
    y = np.asarray(y)
    
    #guessing starting point
    #calculating actual gaussian
    popt,pcov = curve_fit(double_Gauss, x, y, p0=start)
    #recieveing Errors
    #print popt, pcov
    try:
        error=np.sqrt(np.diag(pcov))
    except:
        error=[]
        
    #returns array containing two arrays 1. values 2. corresponding error ( [ [ a , mu , sigma , y0 ] , [ e_a, e_mu , e_ sigma , e_y0 ] ] )
    return [popt,error]
    
    
    
    
    
    
    
    
def Fit_Gaussian(x,y,background=True):
    
    #C reates gaussian Fit
    x = np.asarray(x)
    y = np.asarray(y)
    
    #guessing starting point
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
    
    #calculating actual gaussian
    if background:
        popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma,0])
    else:
        popt,pcov = curve_fit(Gauss_nobackground, x, y, p0=[max(y), mean, sigma])
    #recieveing Errors
    error=np.sqrt(np.diag(pcov))
    
    if not background:
        popt =np.asarray([popt[0],popt[1],popt[2],0])
        error =np.asarray([error[0],error[1],error[2],0])
        
        
        
    #returns array containing two arrays 1. values 2. corresponding error ( [ [ a , mu , sigma , y0 ] , [ e_a, e_mu , e_ sigma , e_y0 ] ] )
    return [popt,error]
    
    
    
    
def Test_gaus_plot(x,y):
    ax=plt.figure("raw values")
    plt.plot(x,y,color="b")
    values , e_values = Fit_Gaussian(x,y,False)
    values2 , e_values2 = Fit_Gaussian(x,y,True)
    print values , e_values
    print Get_pik_counts(x,y,False)
    print values2 , e_values2
    print Get_pik_counts(x,y,True)
        
    y2=[]
    y3=[]
    ax=plt.figure("Gaussian_Fit")
    for channel in x:
        y2.append(Gauss_nobackground(channel, values[0], values[1], values[2]))
        y3.append(Gauss(channel, values2[0], values2[1], values2[2] , values2[3]))
    plt.plot(x,y,color="b")
    plt.plot(x,y2,color="r")
    plt.plot(x,y3,color="g")
    plt.show()
    
def Test_double_gaus_plot(x,y,mean1,mean2,a1,a2):
    ax=plt.figure("raw values")
    plt.plot(x,y,color="b")
    values , e_values = Fit_double_Gaussian(x,y,mean1,mean2,a1,a2)
    print values[0],values[1],values[2]
    print values[3],values[4],values[5]
    print values[6] , values[7]
    y2=[]
    x2=[]
    for i in range(101):
        x2.append(x[0]+(x[-1]-x[0])/100.0*i)
    ax=plt.figure("Gaussian_Fit")
    for channel in x2:
        y2.append(double_Gauss(channel, values[0], values[1], values[2], values[3], values[4], values[5] , values[6], values[7]))
        #y2.append(double_Gauss(channel, 10000, 83, 0.20, 10000, 87, 0.20 , 2000))
    plt.plot(x,y,color="b")
    plt.plot(x2,y2,color="r")
    plt.show()
    
    
    

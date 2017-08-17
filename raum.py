# -*- coding: utf-8 -*-
"""
Created on Tue Jul 04 11:05:15 2017

@author: David
"""
import math
import numpy as np
import time
import matplotlib.pyplot as plt
from pylab import *

def getvariables(x,y,z):
    r=math.sqrt(x**2+y**2+z**2)
    r2=math.sqrt(y**2+z**2)
    phi=math.atan(r2/float(r))
    return [r,phi]
    
def getbeta(x,phi,alpha,r2):
    try:
        ha=math.sqrt(r2**2-(x*math.tan(abs(alpha)))**2)*math.cos(phi)
        return math.atan(ha/math.sqrt(x**2+x**2*math.tan(abs(alpha))**2))
    except :
        print x,phi,alpha,r2,math.atan(r2/float(x))
        math.sqrt(-1)
def getflaeche(beta,alpha):
    #flaeche=alpha*beta*(1+math.sin(beta))
    #flaeche=abs(alpha*beta*2)
    step=0.001
    c=0
    summ=0
    while c<beta:
        summ+=math.acos(math.cos(alpha)*math.cos(c)**2+math.sin(c)**2)*step
        c+=step
    return summ*4
def getint(alpha,start,stop):
    step=(stop-start)/float(int((stop-start)/0.001))
    c=start
    summ=0
    while c<stop:
        summ+=math.acos(math.cos(alpha)*math.cos(c)**2+math.sin(c)**2)*step
        c+=step
    return summ
    


def getmaxalpha(x,r2):
    return abs(math.atan(r2/float(abs(x))))
    
    
def create_template(err,teta=[0,math.pi/2.0],r=[3,5],x=[0.1,15]):
    output=open("raumtemplate.txt","w")
    lines=[]
    t=time.time()
    for i in range(int((teta[1]-teta[0])/math.pi*180)):
        
            lines.append(str(i+teta[0])+";;")
            print i , time.time()-t
            for l in range(int(10*(x[1]-x[0]))):
                
                    lines[-1]+=str(l/10.0+x[0])+",,"
                    for k in range(int(10*(r[1]-r[0]))):
                        alphamax=getmaxalpha(l/10.0+x[0],k/10.0+r[0])
                        flaeche=0
                        error=0
                        for alpha1 in range(100):
                            alp=0+alphamax/100.0*(alpha1+0.5)
                            
                            if abs(alp)>abs(alphamax):
                                print "error",alp, alphamax
                            beta=getbeta(abs(l/10.0+x[0]),(i+teta[0])*math.pi/180.0,float(alp),k/10.0+r[0])
                            flaech=getflaeche(beta,float(alphamax)/100.0)
                            flaeche+=flaech
                            error+=(err[int(beta/math.pi*180)]*flaech)**2
                        
                        lines[-1]+=str(k/10.0+r[0])+":"+str(flaeche*2)+":"+str(math.sqrt(error))+","
                    lines[-1]+=";"
            lines[-1]+="\n"
    output.writelines(lines)
    output.close()
                    
    
    
def getraumwinkel(x,y,z,r2,err):
    r,phi= getvariables(x,y,z)
    alphamax=getmaxalpha(x,r2)
    flaeche=0
    error=0
    for alpha1 in range(100):
        alp=0+alphamax/100.0*(alpha1+0.5)
        
        if abs(alp)>abs(alphamax):
            print "error",alp, alphamax
        beta=getbeta(abs(x),phi,float(alp),r2)
        flaech=getflaeche(beta,float(alphamax)/100.0)
        flaeche+=flaech
        error+=err[int(beta/math.pi*180)]*flaech
    return [2*flaeche,2*error]
    
    
def getraumwinkel3(x,y,z,r2,err,temp):
    r,phi= getvariables(x,y,z)
    return temp.getvalue(phi,r,r2)
    
    
    
def getraumwinkel2(x,y,z,r2):
    r,phi= getvariables(x,y,z)
    alphamax=getmaxalpha(x,r2)
    flaeche=0
    for alpha1 in range(100):
        alp=0+alphamax/100.0*(alpha1+0.5)
        
        if abs(alp)>abs(alphamax):
            print "error",alp, alphamax
        beta=getbeta(abs(x),phi,float(alp),r2)
        flaech=getflaeche(beta,float(alphamax)/100.0)
        flaeche+=flaech
    return flaeche
    
def erfuellt(x,y,z,erfac=[9.0,6.0,6.0],x0=-10.0):
    #drehmatrix=np.asarray([[math.cos(phi),math.sin(phi),0],[-math.sin(phi),math.cos(phi),0],[0,0,1]])
   # x2,y2,z2=np.dot(drehmatrix,np.asarray([x,y,z]))
    a=x>=x0
    b=math.sqrt(((x-x0)/erfac[0])**2+(y/erfac[1])**2+(z/erfac[2])**2)<=1
    return a and b
    
    
    
    
    
def gettotalraumwinkel(x0,x1,y0,y1,z0,z1,r2,n1=105,n2=130,n3=130):
    t=time.time()
    winkel=[]
    current=0
    phi=math.tan(30/180.0*math.pi)
    coss=math.cos(phi)
    sinn=math.sin(phi)
    errors=geterrorarray()
    errorsumme=[]
    c=0
    for x2 in range(n1):
        x=x0+(x1-x0)/float(n1)*(x2+0.5)
        if int(100*c/float(n1))>current:
            current=int(100*c/float(n1))
            print current, "% done ",time.time()-t,"seconds since start , estimated",(time.time()-t)*(100/float(current)-1),"seconds till end current ",len(winkel),"points found"
        c+=1
        for y2 in range(n2):
            y=y0+(y1-y0)/float(n2)*(y2+0.5)
            for z2 in range(n3):
                    z=z0+(z1-z0)/float(n3)*(z2+0.5)
                    if erfuellt(x,y,z):
                        
                        
                        wink,err=getraumwinkel(x*coss+y*sinn,-sinn*x+coss*y,z,r2,errors)
                        errorsumme.append(float(err))
                        winkel.append(wink)
                        
                        
                        
    fehler=0
    for i in range(len(winkel)):
        fehler+=(errorsumme[i]*winkel[i]/float(len(winkel)))**2
    fehler=math.sqrt(fehler)
    print np.mean(winkel),fehler
    return np.mean(winkel)
    
def gettotalraumwinkel2(x0,x1,y0,y1,z0,z1,r2,n1=105,n2=130,n3=130,rfac=[9.0,6.0,6.0]):
    t=time.time()
    winkel=[]
    error=[]
    current=0
    phi=math.tan(30/180.0*math.pi)
    coss=math.cos(phi)
    sinn=math.sin(phi)
    errors=geterrorarray()
    c=0
    aa=template()
    for x2 in range(n1):
        x=x0+(x1-x0)/float(n1)*(x2+0.5)
        if int(100*c/float(n1))>current:
            current=int(100*c/float(n1))
            print current, "% done ",time.time()-t,"seconds since start , estimated",(time.time()-t)*(100/float(current)-1),"seconds till end current ",len(winkel),"points found"
        c+=1
        for y2 in range(n2):
            y=y0+(y1-y0)/float(n2)*(y2+0.5)
            for z2 in range(n3):
                    z=z0+(z1-z0)/float(n3)*(z2+0.5)
                    if erfuellt(x,y,z,rfac):
                        
                        
                        wink,err=getraumwinkel3(x*coss+y*sinn,-sinn*x+coss*y,z,r2,errors,aa)
                        winkel.append(wink)
                        error.append(err)
                        
                        
    print np.mean(winkel),np.mean(error)/math.sqrt(len(error)-1)
    return [np.mean(np.asarray(winkel)),np.mean(error)/math.sqrt(len(error)-1)]
"""
alphamax=200
for alpha1 in np.arange(0+alphamax/200.0,alphamax*1.005,alphamax/100.0):
    print alpha1
print math.atan(3.5/10.0)*180/math.pi
print getraumwinkel(1,0,0,math.tan(30*math.pi/180.0))    
print np.arange(0,1,0.1)
alphamax=10
for alpha1 in np.arange(0+alphamax/200.0,alphamax*1.005,alphamax/100.0):
        print alpha1
        """
        
        
def geterrorwinkel(phi):
    a=4*math.pi*(math.sin(phi/2.0)**2)
    b=getraumwinkel2(1,0,0,math.tan(phi))
    if a!=0:
        return (a-b)/a
    else:
        return 0
def geterrorarray():
    dic={}
    for i in range(91):
        dic.update({i:geterrorwinkel(i/180.0*math.pi)})
    return dic
        
class template():
    def __init__(self):
        temp=[]
        data=open("raumtemplate.txt","r")
        lines=data.readlines()
        for line in lines:
            c=float(line.split(";;")[0])
            temp.append([c])
            for phi in line.split(";;")[1].split(",;"):
                if phi.find("\n")==-1:
                    d=phi.split(",,")[0]
                    temp[-1].append([float(d)])
                    for r in phi.split(",,")[1].split(","):
                        temp[-1][-1].append([float(r.split(":")[0]),float(r.split(":")[1]),float(r.split(":")[2])])
        self.temp=temp
        for i in temp:
            print i[0]
        
    def getstart(self,phi,r,r2):
        i=0
        while self.temp[i][0]<int(phi+0.5):
            i+=1
        
        j=1
        while self.temp[i][j][0]<int(r/0.2+0.5)*0.2:
            j+=1
        
        k=1
        while self.temp[i][j][k][0]<int(r2/0.2+0.5)*0.2:
            k+=1
        return [i,j,k]
        
        
    def getvalue(self,phi,r,r2):
        phi0, r0,r20 = self.getstart(phi*180/math.pi,r,r2)
        dphi=phi-self.temp[phi0][0]/180.0*math.pi
        dr=r-self.temp[phi0][r0][0]
        dr2=r2-self.temp[phi0][r0][r20][0]
        phistep=1
        rstep=1
        r2step=1
        f0=self.temp[phi0][r0][r20][1]
        abl=(self.temp[phi0+phistep][r0][r20][1]*dphi+self.temp[phi0][r0+rstep][r20][1]*dr+self.temp[phi0][r0][r20+r2step][1]*dr2)/2.0
        abl+=(-self.temp[phi0-phistep][r0][r20][1]*dphi-self.temp[phi0][r0-rstep][r20][1]*dr-self.temp[phi0][r0][r20-r2step][1]*dr2)/2.0




        abl2=((self.temp[phi0+2*phistep][r0][r20][1]+self.temp[phi0-2*phistep][r0][r20][1])*dphi*dphi+(self.temp[phi0][r0+2*rstep][r20][1]+self.temp[phi0][r0-2*rstep][r20][1])*dr*dr+(self.temp[phi0][r0][r20+2*r2step][1]+self.temp[phi0][r0][r20-2*r2step][1])*dr2*dr2)/8.0
        abl2+=(-(self.temp[phi0][r0][r20][1])*(dphi*dphi+dr*dr+dr2*dr2))/4.0
        abl2+=(self.temp[phi0+phistep][r0+rstep][r20][1]*dphi*dr+self.temp[phi0+phistep][r0][r20+r2step][1]*dphi*dr2+dr*dr2*self.temp[phi0][r0+rstep][r20+r2step][1])/4.0
        abl2+=(self.temp[phi0-phistep][r0-rstep][r20][1]*dphi*dr+self.temp[phi0-phistep][r0][r20-r2step][1]*dphi*dr2+dr*dr2*self.temp[phi0][r0-rstep][r20-r2step][1])/4.0
        abl2+=-(self.temp[phi0-phistep][r0+rstep][r20][1]*dphi*dr+self.temp[phi0-phistep][r0][r20+r2step][1]*dphi*dr2+dr*dr2*self.temp[phi0][r0-rstep][r20+r2step][1])/4.0
        abl2+=-(self.temp[phi0+phistep][r0-rstep][r20][1]*dphi*dr+self.temp[phi0+phistep][r0][r20-r2step][1]*dphi*dr2+dr*dr2*self.temp[phi0][r0+rstep][r20-r2step][1])/4.0

        ef0=self.temp[phi0][r0][r20][2]**2
        
        
        
        eabl=((self.temp[phi0+phistep][r0][r20][2]*dphi)**2+(self.temp[phi0][r0+rstep][r20][2]*dr)**2+(self.temp[phi0][r0][r20+r2step][1]*dr2)**2)/4.0
        eabl+=((-self.temp[phi0-phistep][r0][r20][2]*dphi)**2+(-self.temp[phi0][r0-rstep][r20][2]*dr)**2+(-self.temp[phi0][r0][r20-r2step][2]*dr2)**2)/4.0
        
        
        eabl2=((self.temp[phi0+2*phistep][r0][r20][2]**2+self.temp[phi0-2*phistep][r0][r20][2]**2)*dphi*dphi*dphi*dphi+(self.temp[phi0][r0+2*rstep][r20][2]**2+self.temp[phi0][r0-2*rstep][r20][2]**2)*dr*dr*dr*dr+(self.temp[phi0][r0][r20+2*r2step][2]**2+self.temp[phi0][r0][r20-2*r2step][2]**2)*dr2*dr2*dr2*dr2)/64.0
        eabl2+=((-(self.temp[phi0][r0][r20][2]))**2*(dphi*dphi*dphi*dphi+dr*dr*dr*dr+dr2*dr2*dr2*dr2))/16.0
        
        
        eabl2+=((self.temp[phi0+phistep][r0+rstep][r20][2]*dphi*dr)**2+(self.temp[phi0+phistep][r0][r20+r2step][2]*dphi*dr2)**2+(dr*dr2*self.temp[phi0][r0+rstep][r20+r2step][2])**2)/16.0
        eabl2+=((self.temp[phi0-phistep][r0-rstep][r20][2]*dphi*dr)**2+(self.temp[phi0-phistep][r0][r20-r2step][2]*dphi*dr2)**2+(dr*dr2*self.temp[phi0][r0-rstep][r20-r2step][2])**2)/16.0
        eabl2+=((self.temp[phi0-phistep][r0+rstep][r20][2]*dphi*dr)**2+(self.temp[phi0-phistep][r0][r20+r2step][2]*dphi*dr2)**2+(dr*dr2*self.temp[phi0][r0-rstep][r20+r2step][2])**2)/16.0
        eabl2+=((self.temp[phi0+phistep][r0-rstep][r20][2]*dphi*dr)**2+(self.temp[phi0+phistep][r0][r20-r2step][2]*dphi*dr2)**2+(dr*dr2*self.temp[phi0][r0+rstep][r20-r2step][2])**2)/16.0
        return [f0 + abl + abl2,math.sqrt(ef0+eabl+eabl2)]
        
        
"""      
x=[]
y=[]
for i in range(180):
    x.append(i/180.0*math.pi)
    y.append(geterrorwinkel(i/180.0*math.pi))
plt.plot(x,y)
"""

def do_everything(template=True):
    if template:
        create_template(geterrorarray())
    a= gettotalraumwinkel2(-10.5,0,-6.5,6.5,-6.5,6.5,3.5,rfac=[9.,6.0,6.0])
    b=gettotalraumwinkel2(-10.5,0,-6.5,6.5,-6.5,6.5,3.5,rfac=[9.1,6.0,6.0])
    c=gettotalraumwinkel2(-10.5,0,-6.5,6.5,-6.5,6.5,3.5,rfac=[8.9,6.0,6.0])
    d=gettotalraumwinkel2(-10.5,0,-6.5,6.5,-6.5,6.5,3.5,rfac=[9.,6.1,6.0])
    e=gettotalraumwinkel2(-10.5,0,-6.5,6.5,-6.5,6.5,3.5,rfac=[9.,5.9,6.0])
    f=gettotalraumwinkel2(-10.5,0,-6.5,6.5,-6.5,6.5,3.5,rfac=[9.,6.0,6.1])
    g=gettotalraumwinkel2(-10.5,0,-6.5,6.5,-6.5,6.5,3.5,rfac=[9.,6.0,5.9])
    
    erros=[(b[0]-c[0])/0.2,(d[0]-e[0])/0.2,(f[0]-g[0])/0.2]
    datei=open("Raumwinkel.txt","w")
    lines=[str(a[0])+"\n",str(np.sqrt(a[1]**2+(0.2*erros[0])**2+(0.2*erros[1])**2+(0.2*erros[2])**2))]
    datei.writelines(lines)
    datei.close()
    print a
#a=template()
#a.getvalue(4,5.0,3.6)
        
        
        
        
        
        
        
        
        
        
        
        
        
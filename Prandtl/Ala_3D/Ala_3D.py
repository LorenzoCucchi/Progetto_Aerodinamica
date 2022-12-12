import numpy as np
import matplotlib.pyplot as plt

#import user functions for the calculations
from single_attitude import single_attitude

class Ellisse:
    def __init__(self, Uinf, rho, bwing, a, b, S, AR, LAR, c1, c, cl_a, alfa_g, alfa0):
        self.Uinf = Uinf
        self.rho=rho 
        self.bwing=bwing
        self.a=a
        self.b=b
        self.S=S
        self.AR=AR
        self.LAR=LAR
        self.c1=c1
        self.c=c


        
#GENERAL SETTINGS

#plots parameters
N_plot=500
Ellisse_calc=True

if Ellisse_calc == True:
    
    # ALA ELLITTICA SETTINGS
    #general settings
    geom_plots=False
    distribution_plots=True
    
    #field parameters
    Ellisse.Uinf=51
    Ellisse.rho=1.225
    
    #Wing parameters
    Ellisse.bwing=11
    Ellisse.a=Ellisse.bwing/2
    Ellisse.b=1
    
    Ellisse.S=np.pi * Ellisse.a*Ellisse.b
    Ellisse.AR=Ellisse.b**2 / Ellisse.S
    Ellisse.LAR=2
    Ellisse.c1=1
    def c(z):
        v=np.zeros(len(z))
        for i in range(len(z)):
            v[i] = Ellisse.c1 * (1 + Ellisse.LAR)* np.sqrt(Ellisse.b**2 - (z[i]**2 / Ellisse.a**2 * Ellisse.b**2))
        return v
    Ellisse.c=c

    Ellisse.cl_a=6
    def alfa_g(z):
        v=np.zeros(len(z))
        for i in range(len(z)):
            v[i] = np.deg2rad(2) 
        return v

    def alfa_0(z):
        v=np.zeros(len(z))
        for i in range(len(z)):
            v[i] = np.deg2rad(-1)
        return v
    
    Ellisse.alfa_g=alfa_g
    Ellisse.alfa_0=alfa_0

    #Geometry plots
    if geom_plots==True:
        zplot=np.linspace(-Ellisse.bwing/2,Ellisse.bwing/2,1000)

        fig, geom = plt.subplots()
        geom.plot(zplot, Ellisse.c(zplot)/(2*Ellisse.c([0])),c='r')
        geom.plot(zplot, -Ellisse.LAR*Ellisse.c(zplot)/(2*Ellisse.c([0])),c='r')
        geom.axis('equal')
        plt.show()

    #Single attitude coefficients and forces calculation
    single_attitude(Ellisse, N_plot, distribution_plots)

    #coefficients with alpha
    

    
   
    


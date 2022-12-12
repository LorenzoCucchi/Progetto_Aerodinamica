import numpy as np
import matplotlib.pyplot as plt

#Calculates the Lift, CL, CDi and Di of the 3D wing at a single alpha

def single_attitude(Ellisse, N_plot, flag_plot):
    from Chev_coeffs import Chev_coeffs #import the script that calculates the values of the Chebycev coefficients
    #mapping
    def z(th):
        v=np.zeros(len(th))
        for i in range(len(th)):
            v[i] = -Ellisse.bwing/2 * np.cos(th)
        return v
    
    N=100 #degree of the Chebycev interpolation
    B = Chev_coeffs(Ellisse,z,N)

    z_v=np.linspace(-Ellisse.bwing/2,Ellisse.bwing/2,N_plot)

    Gamma=np.zeros(N_plot)

    for k in range(N_plot):
        summation=0
        for n in range(N):
            nfor=n+1
            summation=summation + (B[n]*np.sin(nfor*np.arccos(-2/Ellisse.bwing * z_v[k])))
        Gamma[k]=(2*Ellisse.bwing* Ellisse.Uinf* summation)

    #Lift
    CL=-np.pi* Ellisse.AR * B[0]
    L=0.5* Ellisse.rho* Ellisse.Uinf**2 * Ellisse.S * CL

    #induced drag
    delta=0
    for n in range(2,N+1):
        delta=delta + n*(B[n-1]/B[0])**2

    CDi=(CL**2 /(np.pi * Ellisse.AR)) * (1+delta)
    Di= 0.5* Ellisse.rho* Ellisse.Uinf**2 * Ellisse.S * CDi

    if flag_plot==True:

        #2D plot of the distribution of lift
        l=-Gamma*Ellisse.rho*Ellisse.Uinf
        fig, distribution = plt.subplots()
        distribution.plot(z_v, l)
        distribution.set_title('Lift distribution over the wing')
        plt.show()

        zplot=np.linspace(-Ellisse.bwing/2,Ellisse.bwing/2,1000)



        #3D plot of the distribution of lift over the wing plan for an elliptic wing
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        
        #Geometry in the 3D space
        xline=zplot
        yline1= Ellisse.c(zplot)/(2*Ellisse.c([0]))
        zline1=yline1*np.tan(Ellisse.alfa_g([z_v])) #adding the angle of attack
        ax.plot3D(xline,yline1,zline1,'red')
        
        yline2= -Ellisse.LAR*Ellisse.c(zplot)/(2*Ellisse.c([0]))
        zline2= yline2*np.tan(Ellisse.alfa_g([z_v])) #adding the angle of attack
        ax.plot3D(xline,yline2,zline2,'red')

        #Lift distribution in the 3D space
        xlift=z_v
        ylift=np.zeros(len(z_v))
        zlift=l/(Ellisse.S*10) #scaled to fit the plot  l[int(N_plot/2)]
        ax.plot3D(xlift,ylift,zlift,'blue')
        ax.text(0, 0, l[int(N_plot/2)]/(Ellisse.S*10) +1,'CDi='+str(CDi[0]), size=10, ha='center', va='center')
        ax.text(0, 0, l[int(N_plot/2)]/(Ellisse.S*10) +2,'CL='+str(CL[0]), size=10, ha='center', va='center')
        ax.axis('equal')
        plt.show()

    return 


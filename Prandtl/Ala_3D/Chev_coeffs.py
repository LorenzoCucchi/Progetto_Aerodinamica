import numpy as np

#calculates the chebicev coefficients to numerically calculate Gamma

def Chev_coeffs(Ellisse,z,N):

    b=Ellisse.bwing
    c=Ellisse.c
    cl_a=Ellisse.cl_a
    alfa_g=Ellisse.alfa_g
    alfa_0=Ellisse.alfa_0

    #there is probably a better way to work with indexes starting from 0 but idk
    thi=np.zeros(N)
    for i in range(1,N+1):
        thi[i-1]=np.pi/N * (i-0.5)
    
    A=np.zeros((N,N))
    F=np.zeros((N,1))
    for i in range(1,N+1):
        for n in range(1,N+1):
            A[i-1,n-1]=( -4*b*np.sin(n*thi[i-1]) / c(z([thi[i-1]])) ) - (n*cl_a* np.sin(n*thi[i-1])/np.sin(thi[i-1]))
        
        F[i-1]=cl_a*( alfa_g(z([thi[i-1]])) - alfa_0(z([thi[i-1]])) )

    B=np.linalg.solve(A, F)

    return B
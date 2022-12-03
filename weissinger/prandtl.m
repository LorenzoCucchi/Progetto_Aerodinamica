%% Calcolo CL in un Cesna 172
clear all
close all
%% DATI
b=11; %wingspan [m]
S=16.2; %Wing area [m^2]
lamda=b^2/S; %Aspect Ratio
c=S/b; % corda [m]; ipotizzo per semplicità un area quadrata 
U=226/3.6; %[m/s] velodcità di crociera
h=4100; %altezza di crociera [m]
[T,a,P,rho]=atmosisa(h);
mu=1.81*10^(-6);
Re=rho*U*c/mu;
%nel caso viscoso, il flusso è attaccato si può verificare con Xfoil
W_vuoto=502*9.8; %[kg]
W_pieno=757*9.8; %[kg]
% Angolo di calettamento (Washout)
% NACA 2412
m=0.02;
p=0.4;
alfa_0=Calcolo_Alfa0(m,p);
alfa_0=0;
Cl_a=2*pi;
%% Lifting line Theory
num=14;
agmax=7;
agmin=5;
alfa_g=0*pi/180;
%a_g=@(z) (-4.*(agmax-agmin/2)./b.^2.*z.^2+agmax).*pi/180;
%Corda=@(z) 2.*S./pi./b.*sqrt(1-4.*z.^2./b.^2);
alfa_vect=linspace(-3,10,num);
for t=1:num
alfa=alfa_vect(t)*pi/180;

%Punti di Chebyshev
N=10;
for i=1:N
    teta(i)=pi/N*(i-1/2);
    f(i)=Cl_a.*(alfa_g+alfa-alfa_0);
end
for i=1:N
    for n=1:N
        A(i,n)=(-4*b/c-Cl_a*n/sin(teta(i)))*sin(n*teta(i));
    end
end
B=A\f';

CL(t)=-pi*B(1)*lamda;
L=0.5*rho*U^2*CL*S;
CD(t)=pi*lamda*sum(B.^2);
D=0.5*rho*U^2*CD*S;
figure (1)
hold on
labels=cellstr(num2str(alfa_vect(t)));
plot(CD(t)+0.033,CL(t),'*','color','g');
%text(CD(t)+0.033,CL(t),labels);
e(t)=CL(t)^2/pi/lamda/CD(t);
hold on
end
figure (1)
plot(CD+0.033,CL,'r','linewidth',1);
xlabel("Cd");
ylabel("Cl");
title("Curva Polare Cesna 172R");
grid on
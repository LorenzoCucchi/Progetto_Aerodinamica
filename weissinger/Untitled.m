clear all
close all
tic
%%
% DATI CESNNA 172
b=11;
Delta=0*pi/180;
d=1*pi/180;
S=16.2;
C_t=1.4;
C_r=S*2/b-C_t;
alfa=6;
beta=0;
N=5;
[~,~,~,rho]=atmosisa(4100);
U_inf=226/3.6;
c_med=(C_r+C_t)/2;
i_w=4;
%% Costruisco la geometria
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0],i_w);
U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
%% Costruisco la geometria dellla coda
% dati cessna 172
S_tail=2;
b_tail=3.4;
C_t_tail=0.4; %DImensione alle estremità
Delta_tail=6/180*pi; % angolo di freccia
d_tail= 0/180*pi ; %Angolo di diedro
C_r_tail=S_tail*2/b_tail-C_t_tail;
c_med_tail=C_r_tail+C_t_tail;
distance=4+C_r;
hh=-0.2;
i_t=1;
% dati pamadi
[X_tail,Y_tail,Z_tail,p_tail,f_tail,M_tail]=Geometria(b_tail,Delta_tail,C_r_tail,C_t_tail,d_tail,N,[distance 0 hh],i_t);
%% calcolo delle due polari e a_t, a_w
%% stack geometrie
p_tot=p;
for i=1:f_tail
    p_tot.panels(i+f)=p_tail.panels(i);
end
f_tot=f_tail+f;

%% Calcolo distribuzione di vorticità
[Gamma,A1,b1,FX,FY,FZ]=LinearSystem(p,f,U);
[Gamma_alone_tail,A3,b3,FX3,FY3,FZ3]=LinearSystem(p_tail,f_tail,U);
[Gamma2,A2,b2,FX2,FY2,FZ2]=LinearSystem(p_tot,f_tot,U);
Q=11;
alfa_vect=linspace(-5,5,Q);
for t=1:Q
    b_wing_tail=zeros(f_tot,1);
    alfa=alfa_vect(t);
    U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
   
    % calcolo circolazione totale
    for i=1:f_tot
        Normal=p_tot.panels(i).n;
        b_wing_tail(i)=-dot(U,Normal);
    end
    Gamma2=A2\b_wing_tail;
    Gamma_wing_d=Gamma2(1:2*M*N);
Gamma_tail_d=Gamma2(2*M*N+1:end);
%%
Gamma_matrix_d=[];
% calcolo per quelle accoppiate
for k=0:2*N-1
    Gamma_matrix_d=[Gamma_matrix_d;Gamma_tail_d((1+k*M_tail):(M_tail+k*M_tail))'];
end
Gamma_matrix_d=[fliplr(Gamma_matrix_d(N+1:end,:)) Gamma_matrix_d(1:N,:)];
Gamma_matrix1_d=[];
for k=0:2*N-1
    Gamma_matrix1_d=[Gamma_matrix1_d;Gamma_wing_d((1+k*M):(M+k*M))'];
end
Gamma_matrix1_d=[fliplr(Gamma_matrix1_d(N+1:end,:)) Gamma_matrix1_d(1:N,:)];
[Cl_2D_s,L_tail_d]=portanza(M_tail,Gamma_matrix_d,rho,U_inf,c_med_tail,d_tail,b_tail,N);
[Cl_2D_d,L_wing_d]=portanza(M,Gamma_matrix1_d,rho,U_inf,c_med,d,b,N);
Cl_tot(t)=(L_tail_d+L_wing_d)/(0.5*rho*U_inf^2*S);
CL(t)=0.0805*(alfa+i_w)+S_tail/S*0.0761*(alfa*(1-0.3)+i_t-0.5);
end
figure ()
plot(alfa_vect,Cl_tot,'b','linewidth',1);
hold on
plot(alfa_vect,CL,'r','linewidth',1);
grid on
legend("CL wessinger","CL metodo 2 superfici");
title('Modello due superfici')
xlabel("\alpha");
ylabel ("CL")
clear all
tic
%% 
%% Input iniziali
C_r =2; %Dimensione corda alla basa
C_t = 1.5; %DImensione alle estremità
b = 10;       % b
Delta =20/180*pi; % angolo di freccia
d = 5/180*pi ; %Angolo di diedro
N=10; %Numero di pannelli

C_r =2; %Dimensione corda alla basa
C_t = 1.5; %DImensione alle estremità
b = 10;       % b
Delta =20/180*pi; % angolo di freccia
d = 5/180*pi ; %Angolo di diedro

% DATI B737
C_r=18*0.3048+8.75*0.0254; %Dimensione corda alla basa
C_t=4*0.3048+1.25*0.0254; %DImensione alle estremità
b =112*0.3048+7*0.0254;       % b
Delta =25/180*pi; % angolo di freccia
d = 6/180*pi ; %Angolo di diedro
%N=5; %Numero di pannelli

N=5; %Numero di pannelli
alfa=4;
beta=5;
rho=1;
U_inf=1;
c_med=(C_r+C_t)/2;
S=b*(C_r+C_t)/2;
%% Costruisco la geometria
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0],0);
U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
%% Costruisco la geometria dellla coda
% DATI B737
S_tail=0.092903*352.8;
C_t_tail=2*0.3048+6.6*0.0254; %DImensione alle estremità
b_tail =47*0.3048+1*0.0254;       % b
Delta_tail=30/180*pi; % angolo di freccia
d_tail= 7/180*pi ; %Angolo di diedro
C_r_tail=S_tail*2/b_tail-C_t_tail;
Cale=-4; %Calettamento
distance=C_r+b/3;
G=[20 0 0];
[X_tail,Y_tail,Z_tail,p_tail,f_tail,M_tail]=Geometria(b_tail,Delta_tail,C_r_tail,C_t_tail,d_tail,N,[distance 0 1],Cale);
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
%% rappresentazione grafica del risultato della funzione
Gamma_wing=Gamma2(1:2*M*N);
Gamma_tail=Gamma2(2*M*N+1:end);
%% riscrivo gamma come una matrice n*2M
%% riscrivo gamma come una matrice n*2M
Gamma_matrix_wing=[];
for f=0:2*N-1
    Gamma_matrix_wing=[Gamma_matrix_wing;Gamma((1+f*M):(M+f*M))'];
end
Gamma_matrix_wing=[fliplr(Gamma_matrix_wing(N+1:end,:)) Gamma_matrix_wing(1:N,:)];
Gamma_matrix_wing_d=[];
for f=0:2*N-1
    Gamma_matrix_wing_d=[Gamma_matrix_wing_d;Gamma2((1+f*M):(M+f*M))'];
end
Gamma_matrix_wing_d=[fliplr(Gamma_matrix_wing_d(N+1:end,:)) Gamma_matrix_wing_d(1:N,:)];
Gamma_m_alone_tail=[];
for f=0:2*N-1
    Gamma_m_alone_tail=[Gamma_m_alone_tail;Gamma_alone_tail((1+f*M_tail):(M_tail+f*M_tail))'];
end
Gamma_m_alone_tail=[fliplr(Gamma_m_alone_tail(N+1:end,:)) Gamma_m_alone_tail(1:N,:)];
Gamma_m_tail=[];
for f=0:2*N-1
    Gamma_m_tail=[Gamma_m_tail;Gamma_tail((1+f*M_tail):(M_tail+f*M_tail))'];
end
Gamma_m_tail=[fliplr(Gamma_m_tail(N+1:end,:)) Gamma_m_tail(1:N,:)];
%% Calcolo portanza 
[Cl_2D_s,L_s]=portanza(M,Gamma_matrix_wing,rho,U_inf,c_med,d,b,N);
[Cl_2D_d,L_d]=portanza(M,Gamma_matrix_wing_d,rho,U_inf,c_med,d,b,N);
Cl_S_wing=L_s/(0.5*rho*U_inf^2*S)
Cl_d_wing=L_d/(0.5*rho*U_inf^2*S)
[Cl_2D_s_tail,L_s_tail]=portanza(M_tail,Gamma_m_alone_tail,rho,U_inf,c_med,d_tail,b_tail,N);
[Cl_2D_d_tail,L_d_tail]=portanza(M_tail,Gamma_m_tail,rho,U_inf,c_med,d_tail,b_tail,N);
Cl_S_tail=L_s_tail/(0.5*rho*U_inf^2*S_tail)
Cl_d_tail=L_d_tail/(0.5*rho*U_inf^2*S_tail)

y_vect=linspace(-b/2,b/2,2*M);
figure()
plot(y_vect,Cl_2D_s,'b','linewidth',1);
hold on
plot(y_vect,Cl_2D_d,'r','linewidth',1);
legend("Ala singola","Ala con coda");
grid on
title (['CL_2_d con ', '\alpha',' = ', num2str(alfa),'° ',...
    '\beta',' = ', num2str(beta),'°',' Calettamento',' = ', num2str(Cale),'°',' Distanza',' = ', num2str(distance),'m'])
y_vect_tail=linspace(-b_tail/2,b_tail/2,2*M_tail);
figure()
plot(y_vect_tail,Cl_2D_s_tail,'b','linewidth',1);
hold on
plot(y_vect_tail,Cl_2D_d_tail,'r','linewidth',1);
legend("Coda singola","Coda con ala");
grid on
title (['CL_2_d con ', '\alpha',' = ', num2str(alfa),'° ',...
    '\beta',' = ', num2str(beta),'°',' Calettamento',' = ', num2str(Cale),'°',' Distanza',' = ', num2str(distance),'m'])

%% Calcolo velocità indotta
Time=toc;
fprintf("Time: %.3f s \n",Time);

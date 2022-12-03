clear all
close all
tic

%% dati UAV
% DATI CESNNA 172
b=11;
Delta=0*pi/180;
d=1*pi/180;
S=16.2;
C_t=1.4;
C_r=S*2/b-C_t;
alfa=0;
beta=0;
rho=1.225;
U_inf=1;
c_med=(C_r+C_t)/2;
i_w=4;
N=5;
%% Costruisco la geometria
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0],i_w);
U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
%% Costruisco la geometria dellla coda
%% 
C_r_tail=0.65; %DImensione alle estremità
C_t_tail=0.65;
b_tail=2.5;       % b
S_tail=1.625;
Delta_tail=0; % angolo di freccia
d_tail= 0; %Angolo di diedro
c_med_tail=(C_t_tail+C_r_tail)/2;
Cale=1;
distance=3.3;
S_tail=2;
b_tail=3.4;
C_t_tail=0.4; %DImensione alle estremità
Delta_tail=6/180*pi; % angolo di freccia
d_tail= 0/180*pi ; %Angolo di diedro
C_r_tail=S_tail*2/b_tail-C_t_tail;
distance=4+C_r;
hh=-0.2;
G=[20 0 0];
[X_tail,Y_tail,Z_tail,p_tail,f_tail,M_tail]=Geometria(b_tail,Delta_tail,C_r_tail,C_t_tail,d_tail,N,[distance 0 hh],Cale);
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
%GraficoIniziale(p,Gamma,U,0)
figure (1)
Gamma_wing=Gamma2(1:2*M*N);
Gamma_tail=Gamma2(2*M*N+1:end);
visual_circ(X,Y,Z,Gamma_wing,M,N)
hold on
visual_circ(X_tail,Y_tail,Z_tail,Gamma_tail,M_tail,N)
title (['\Gamma con 2 sup', '\alpha',' = ', num2str(alfa),'° ',...
    '\beta',' = ', num2str(beta),'°'])
figure(2)
visual_circ(X,Y,Z,Gamma,M,N)
title (['\Gamma con ala singola', '\alpha',' = ', num2str(alfa),'° ',...
    '\beta',' = ', num2str(beta),'°'])
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
[cd,gamma_scia]=Resistenza(rho,U_inf,Gamma_matrix_wing,Y,Z,M,p,alfa,beta,S);
%% Calcolo portanza 
[Cl_2D_s,L_s]=portanza(M,Gamma_matrix_wing,rho,U_inf,c_med,d,b,N);
[Cl_2D_d,L_d]=portanza(M,Gamma_matrix_wing_d,rho,U_inf,c_med,d,b,N);
Cl_S_wing=L_s/(0.5*rho*U_inf^2*S)
Cl_d_wing=L_d/(0.5*rho*U_inf^2*S)
[Cl_2D_s_tail,L_s_tail]=portanza(M_tail,Gamma_m_alone_tail,rho,U_inf,c_med,d_tail,b_tail,N);
[Cl_2D_d_tail,L_d_tail]=portanza(M_tail,Gamma_m_tail,rho,U_inf,c_med,d_tail,b_tail,N);
Cl_S_tail=L_s_tail/(0.5*rho*U_inf^2*S_tail)
Cl_d_tail=L_d_tail/(0.5*rho*U_inf^2*S_tail)
%B737 Prova
h=10668;
V=926.10/3.6;
W=50000*9.8;
[~,~,~,rho2]=atmosisa(h);
Cl2=W/(0.5*rho2*V^2*S);
% Andamento del Cl_2D
y_vect=linspace(-b/2,b/2,2*M);
figure()
plot(y_vect,Cl_2D_s,'b','linewidth',1);
hold on
plot(y_vect,Cl_2D_d,'r','linewidth',1);
legend("Ala singola","Ala con coda");
grid on
title (['CL_2_d con ', '\alpha',' = ', num2str(alfa),'° ',...
    '\beta ',' = ', num2str(beta),'°',' Freccia',' = ', num2str(Delta*180/pi),'°'])
y_vect_tail=linspace(-b_tail/2,b_tail/2,2*M_tail);
xlabel("b [m]");
ylabel(" Cl");
figure()
plot(y_vect_tail,Cl_2D_s_tail,'b','linewidth',1);
hold on
plot(y_vect_tail,Cl_2D_d_tail,'r','linewidth',1);
legend("Coda singola","Coda con ala");
grid on
title (['CL_2_d con ', '\alpha',' = ', num2str(alfa),'° ',...
    '\beta',' = ', num2str(beta),'°',' Calettamento',' = ', num2str(Cale),'°'])
%% Calcolo delle forze e momenti:
% stavo cercando di utilizzare ciò che è scritto nelle vecchie slide, ma
% nulla da fare, bho ho numerosi dubbi
fprintf("Calcolo f_orze.... \n");
b2w=[cos(alfa*pi/180)*cos(beta*pi/180) -sin(beta*pi/180) sin(alfa*pi/180)*cos(beta*pi/180) ;...
           cos(alfa*pi/180)*sin(beta*pi/180)  cos(beta*pi/180) sin(alfa*pi/180)*sin(beta*pi/180) ;...
           -sin(alfa*pi/180)               0         cos(alfa*pi/180)        ];
[F,Moment]=force(p,Gamma,rho,U,G);
F_wind=b2w*F';
CF_wind=F_wind/(0.5*rho*U_inf^2*S);
CM_wind=[Moment(1)/(0.5*rho*U_inf^2*S*b) Moment(2)/(0.5*rho*U_inf^2*S*c_med) Moment(3)/(0.5*rho*U_inf^2*S*b)];
c_med_tail=(C_r_tail+C_t_tail)/2;
[F_tail,Moment_tail]=force(p,Gamma,rho,U,G);
F_wind_tail=b2w*F_tail';
CF_wind_tail=F_wind_tail/(0.5*rho*U_inf^2*S_tail);
CM_wind_tail=[Moment_tail(1)/(0.5*rho*U_inf^2*S_tail*b_tail) Moment_tail(2)/(0.5*rho*U_inf^2*S_tail*c_med_tail) Moment_tail(3)/(0.5*rho*U_inf^2*S_tail*b_tail)];
[F_wing,Moment_wing]=force(p,Gamma,rho,U,G);
F_wind_wing=b2w*F_wing';
CF_wind_wing=F_wind_wing/(0.5*rho*U_inf^2*S);
CM_wind_wing=[Moment_wing(1)/(0.5*rho*U_inf^2*S*b) Moment_wing(2)/(0.5*rho*U_inf^2*S*c_med) Moment_wing(3)/(0.5*rho*U_inf^2*S*b)];
[F_tail_alone,Moment_tail_alone]=force(p,Gamma,rho,U,G);
F_wind_tail_alone=b2w*F_tail_alone';
CF_wind_tail_alone=F_wind_tail_alone/(0.5*rho*U_inf^2*S_tail);
CM_wind_tail_alone=[Moment_tail_alone(1)/(0.5*rho*U_inf^2*S_tail*b_tail) Moment_tail_alone(2)/(0.5*rho*U_inf^2*S_tail*c_med_tail) Moment_tail_alone(3)/(0.5*rho*U_inf^2*S_tail*b_tail)];
% la portanza è quasi dimezzata della coda posteriore
%% Calcolo velocità indotta
Time=toc;
fprintf("Time: %.3f s \n",Time);


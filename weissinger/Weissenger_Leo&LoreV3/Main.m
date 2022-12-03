clear all
close all
tic
%% 
%% Input iniziali
C_r =2; %Dimensione corda alla basa
C_t = 1.5; %DImensione alle estremità
b = 10;       % b
Delta =20/180*pi; % angolo di freccia
d = 5/180*pi ; %Angolo di diedro
N=10; %Numero di pannelli

% DATI B737
C_r=18*0.3048+8.75*0.0254; %Dimensione corda alla basa
C_t=4*0.3048+1.25*0.0254; %DImensione alle estremità
b =112*0.3048+7*0.0254;       % b
Delta =25/180*pi; % angolo di freccia
d = 6/180*pi ; %Angolo di diedro
%N=5; %Numero di pannelli
C_r=0.350; %Dimensione corda alla basa
C_t=0.117; %DImensione alle estremità
b =0.117;       % b
Delta =pi/2-atan(0.117/(0.350-0.117)); % angolo di freccia
d = 0/180*pi ; %Angolo di diedro
%N=5; %Numero di pannelli
N=5; %Numero di pannelli
alfa=5;
beta=0;
rho=1;
U_inf=1;
c_med=(C_r+C_t)/2;
S=b*(C_r+C_t)/2;
%% Costruisco la geometria
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0]);
U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
%% Costruisco la geometria dellla coda
% DATI B737
S_tail=0.092903*352.8;
C_t_tail=2*0.3048+6.6*0.0254; %DImensione alle estremità
b_tail =47*0.3048+1*0.0254;       % b
Delta_tail=30/180*pi; % angolo di freccia
d_tail= 7/180*pi ; %Angolo di diedro
C_r_tail=S_tail*2/b_tail-C_t_tail;
G=[20 0 0];
[X_tail,Y_tail,Z_tail,p_tail,f_tail,M_tail]=Geometria(b_tail,Delta_tail,C_r_tail,C_t_tail,d_tail,N,[C_r+b/3 0 1]);

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
Gamma_matrix=[];
for f=0:2*N-1
    Gamma_matrix=[Gamma_matrix;Gamma((1+f*M):(M+f*M))'];
end
Gamma_matrix=[fliplr(Gamma_matrix(N+1:end,:)) Gamma_matrix(1:N,:)];
Gamma_matrix1=[];
for f=0:2*N-1
    Gamma_matrix1=[Gamma_matrix1;Gamma2((1+f*M):(M+f*M))'];
end
Gamma_matrix1=[fliplr(Gamma_matrix1(N+1:end,:)) Gamma_matrix1(1:N,:)];
%% Calcolo portanza 
[Cl_2D_s,L_s]=portanza(M,Gamma_matrix,rho,U_inf,c_med,d,b,N);
[Cl_2D_d,L_d]=portanza(M,Gamma_matrix1,rho,U_inf,c_med,d,b,N);
Cl=L_s/(0.5*rho*U_inf^2*S);
%B737 Prova
h=10668;
V=926.10/3.6;
W=50000*9.8;
[~,~,~,rho2]=atmosisa(h);
Cl2=W/(0.5*rho2*V^2*S);
% Andamento del Cl_2D
y_vect=linspace(-b/2,b/2,2*M);
figure(3)
plot(y_vect,Cl_2D_s,'b','linewidth',1);
hold on
plot(y_vect,Cl_2D_d,'r','linewidth',1);
legend("Ala singola","Ala con coda");
grid on
title (['CL_2_d con ', '\alpha',' = ', num2str(alfa),'° ',...
    '\beta',' = ', num2str(beta),'°',' Freccia',' = ', num2str(Delta*180/pi),'°'])
%% Calcolo delle forze e momenti:
% stavo cercando di utilizzare ciò che è scritto nelle vecchie slide, ma
% nulla da fare, bho ho numerosi dubbi
fprintf("Calcolo forze.... \n");
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
%% Calcolo polare 
fprintf("Calcolo polare e curva CL_alpha.... \n");
alfa_vect=linspace(-5,10,15);
f=2*M*N;
CF_wind_new=zeros(15,3);
CM_wind_new=zeros(15,3);
for t=1:15
    alfa=alfa_vect(t);
    U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
    b2w=[cos(alfa*pi/180)*cos(beta*pi/180) -sin(beta*pi/180) sin(alfa*pi/180)*cos(beta*pi/180) ;...
           cos(alfa*pi/180)*sin(beta*pi/180)  cos(beta*pi/180) sin(alfa*pi/180)*sin(beta*pi/180) ;...
           -sin(alfa*pi/180)               0         cos(alfa*pi/180)        ];
    b_new=zeros(f,1);
    for i=1:f
        Normal=p.panels(i).n;
        b_new(i)=-dot(U,Normal);
    end
    Gamma_new=A1\b_new;
    [F_new,M_new]=force(p,Gamma_new,rho,U,G);
    F_wind=b2w*F_new';
CF_wind_new(t,:)=F_wind'/(0.5*rho*U_inf^2*S);
CM_wind_new(t,:)=[M_new(1)/(0.5*rho*U_inf^2*S*b) M_new(2)/(0.5*rho*U_inf^2*S*c_med) M_new(3)/(0.5*rho*U_inf^2*S*b)];
end
figure(4) 
plot(alfa_vect,CF_wind_new(:,3),'r','linewidth',1);
title("Curva CL_\alpha");
xlabel("\alpha");
ylabel("CL");
grid on
CL_alpha=(CF_wind_new(8,3)-CF_wind_new(7,3))/(alfa_vect(8)*pi/180-alfa_vect(7)*pi/180);
figure(5)
plot(CF_wind_new(:,1),CF_wind_new(:,3),'r','linewidth',1);
title("Drag-Polar");
xlabel("CF");
ylabel("CL");
grid on

%% Calcolo velocità indotta
Time=toc;
fprintf("Time: %.3f s \n",Time);


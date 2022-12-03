clear all
close all
tic

%% dati UAV
b =7.905;       % b
%S=b^2/7.5;
C_r=1.0543; %Dimensione corda alla basa
C_t=1.0543; %DImensione alle estremità
S=C_r*b
%C_r=S/b; %Dimensione corda alla basa
%C_t=S/b; %DImensione alle estremità
Delta =0/180*pi; % angolo di freccia
d = 4/180*pi ; %Angolo di diedro
N=5; %Numero di pannelli
alfa=4;
beta=0;
rho=1;
U_inf=1;
c_med=(C_r+C_t)/2;
%% Costruisco la geometria
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0],0);
U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
%% Costruisco la geometria dellla coda
%% Calcolo distribuzione di vorticità
[Gamma,A1,b1,FX,FY,FZ]=LinearSystem(p,f,U);
%% rappresentazione grafica del risultato della funzione
%GraficoIniziale(p,Gamma,U,0)
figure(1)
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
%% calcolo resistenza
% HO una matrice Gamma_matrix di dimensione N*2M, calcolo i vortici in scia
% per cui prima sommo tutti i vortici rispetto ad N e ottengo il vortice a
% ferro di cavallo "Totale" visto dal bordo di uscita
Gamma_bordo_uscita=sum(Gamma_matrix_wing,1);
l_2d = rho*U_inf*Gamma_bordo_uscita;
Y_tot=[-fliplr(Y(1,2:end)) Y(1,:)];
Z_tot=[fliplr(Z(1,2:end)) Z(1,:)];
dy = abs(Y_tot(1,1)-Y_tot(1,2));

gamma_scia = [Gamma_bordo_uscita(1),-Gamma_bordo_uscita(1:(end-1))+Gamma_bordo_uscita(2:end),-Gamma_bordo_uscita(end)];
y = Y_tot(1,:);
z = Z_tot(1,:);

y_c = y(1:end-1)+0.5*dy;
z_c = 0.5*(z(1:end-1)+z(2:end));

v_ind = zeros(size(y_c));

for i = 1:length(y_c)
    for t = 1:length(gamma_scia)
        dist = sqrt((y_c(i)-y(t))^2+(z_c(i)-z(t))^2);
        v_ind(i) = v_ind(i) + gamma_scia(t)*(-y_c(i)+y(t))/(4*pi*dist^2);
    end
end
y_vect=zeros(1,M);
for i=1:M
    y_vect(i)=p.panels(i).P(2);
end
y_vect=[-fliplr(y_vect) y_vect];
alfa_ind = -atan(v_ind/U_inf); % vettore incidenze indotte
figure ()
quiver(y_vect,zeros(1,2*M),zeros(1,2*M),v_ind), grid on
title (['Velocita  indotta  ', '\alpha',' = ',...
    num2str(alfa),'°', '\beta',' = ', num2str(beta),'°'])
xlabel 'y', ylabel 'z '
d2d=l_2d.*sin(alfa_ind);
cd=sum(d2d)*dy/(0.5*rho*U_inf^2*S);
gamma_scia = [Gamma_bordo_uscita(1),-Gamma_bordo_uscita(1:(end-1))+Gamma_bordo_uscita(2:end),-Gamma_bordo_uscita(end)];

%% Calcolo delle forze e momenti:
% stavo cercando di utilizzare ciò che è scritto nelle vecchie slide, ma
% nulla da fare, bho ho numerosi dubbi
%% Calcolo polare
fprintf("Calcolo polare e curva CL_alpha.... \n");
alfa_vect=linspace(0,10,11);
f=2*M*N;
CF_wind_new=zeros(11,3);
CM_wind_new=zeros(11,3);
ris=struct;
for t=1:11
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
    [F_new,M_new]=force(p,Gamma_new,rho,U,[0 0 0]);
    F_wind=b2w*F_new';
    
    CF_wind_new(t,:)=F_wind'/(0.5*rho*U_inf^2*S);
    CM_wind_new(t,:)=[M_new(1)/(0.5*rho*U_inf^2*S*b) M_new(2)/(0.5*rho*U_inf^2*S*c_med) M_new(3)/(0.5*rho*U_inf^2*S*b)];
    
    ris.num(t).alfa=alfa;
    ris.num(t).CL=CF_wind_new(t,3);
    ris.num(t).CD=CF_wind_new(t,1);
    ris.num(t).CMa=CM_wind_new(t,2);
end
figure()
plot(alfa_vect,CF_wind_new(:,3),'r','linewidth',1);
title("Curva CL_\alpha");
xlabel("\alpha");
ylabel("CL");
grid on
CL_alpha=(CF_wind_new(8,3)-CF_wind_new(7,3))/(alfa_vect(8)*pi/180-alfa_vect(7)*pi/180);
figure()
plot(CF_wind_new(:,1),CF_wind_new(:,3),'r','linewidth',1);
title("Drag-Polar");
xlabel("CD");
ylabel("CL");
grid on
%% Calcolo velocità indotta
Time=toc;
fprintf("Time: %.3f s \n",Time);
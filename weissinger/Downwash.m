clear all
close all
tic
%%
%% Input iniziali
% DATI CESNNA 172
b=11;
Delta=0*pi/180;
d=1*pi/180;
S=16.2;
C_t=1.4;
C_r=S*2/b-C_t;
N=5; %Numero di pannelli
alfa=0;
beta=0;
rho=1;
U_inf=1;
c_med=(C_r+C_t)/2;
i_w=4;
%% Costruisco la geometria
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0],i_w);
U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)]
a=U(1)
h=U(2);
c=U(3);
U_ind1=[a;h;c-0.0101];
U_ind=norm(U_ind1);
%% Costruisco la geometria dellla coda
% DATI B737

% dati cessna 172
S_tail=2;
b_tail=3.4;
C_t_tail=0.4; %DImensione alle estremità
Delta_tail=6/180*pi; % angolo di freccia
d_tail= 0/180*pi ; %Angolo di diedro
C_r_tail=S_tail*2/b_tail-C_t_tail;
distance=4+C_r;
hh=-0.2;
c_med_tail=(C_r_tail+C_t_tail)/2;
% dati pamadi
AR=b^2/S;
lamda=C_t/C_r;

lh=distance-C_r/4+C_r_tail/4;
K_AR=1/AR-1/(1+AR^1.7);
K_lamda=(10-3*lamda)/7;
K_r=(1-hh/b)/(2*lh/b)^(1/3);
deps=4.44*(K_AR*K_lamda*K_r*sqrt(cos(Delta)))^1.19
i_t=1;

[X_tail,Y_tail,Z_tail,p_tail,f_tail,M_tail]=Geometria(b_tail,Delta_tail,C_r_tail,C_t_tail,d_tail,N,[distance 0 hh],i_t);
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
alfa_vect=linspace(-5,10,30);
for t=1:30
    b_tail_s=zeros(f_tail,1);
    b_tail_wing=zeros(f_tot,1);
    alfa=alfa_vect(t);
    U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];
    % calcolo il vettore termini noti per il sistema composto da coda,
    % dunque ha dimensione N*2*M_tail
    for i=1:f_tail
        Normal=p_tail.panels(i).n;
        b_tail_s(i)=-dot(U,Normal);
    end
    %Calco della cricolazione sulla coda singola
    Gamma_alone_tail=A3\b_tail_s;
    for i=1:f_tot
        Normal=p_tot.panels(i).n;
        b_tail_wing(i)=-dot(U,Normal);
    end
    Gamma_tail_wing=A2\b_tail_wing;
    Gamma_tail_d=Gamma_tail_wing(2*M*N+1:end);
    % Costruisco le matrici
    Gamma_matrix=[];
    for f=0:2*N-1
        Gamma_matrix=[Gamma_matrix;Gamma_alone_tail((1+f*M_tail):(M_tail+f*M_tail))'];
    end
    Gamma_matrix=[fliplr(Gamma_matrix(N+1:end,:)) Gamma_matrix(1:N,:)];
    Gamma_matrix1=[];
    for f=0:2*N-1
        Gamma_matrix1=[Gamma_matrix1;Gamma_tail_d((1+f*M_tail):(M_tail+f*M_tail))'];
    end
    Gamma_matrix1=[fliplr(Gamma_matrix1(N+1:end,:)) Gamma_matrix1(1:N,:)];
    [Cl_2D_s,L_s]=portanza(M_tail,Gamma_matrix,rho,U_inf,c_med_tail,d_tail,b_tail,N);
    [Cl_2D_d,L_d]=portanza(M_tail,Gamma_matrix1,rho,U_inf,c_med_tail,d_tail,b_tail,N);
    Cl_S_alfa(t)=L_s/(0.5*rho*U_inf^2*S_tail);
    Cl_d_alfa(t)=L_d/(0.5*rho*U_inf^2*S_tail);
end
%%
figure ()
plot(alfa_vect,Cl_S_alfa,'b','linewidth',1);
hold on
plot(alfa_vect,Cl_d_alfa,'r','linewidth',1);
grid on
legend("Coda da sola","Coda con ala davanti");
title(['Effetto Downwash:',' \alpha',' = ', num2str(alfa),'°',' i_t',' = ', num2str(i_t),'°',' i_w',' = ', num2str(i_w),'°',' Distanza',' = ', num2str(distance),'m'])
xlabel("\alpha");
%% Studio down
c1=polyfit(alfa_vect,Cl_S_alfa,1);
c2=polyfit(alfa_vect,Cl_d_alfa,1);
m1=c1(1);
m2=c2(1);
q1=c1(2);
q2=c2(2);
Eps_0=(q1-q2)/m2;
Eps_alfa=(m1/m2-1)
Eps=@(a) Eps_0+Eps_alfa.*a;
Cl_s=@(a) m1.*a+q1;
Cl_d=@(a) m2.*a+q2;
errore= (Eps_alfa-deps)/(Eps_alfa+deps)*2*100


%cl_tail_wing(

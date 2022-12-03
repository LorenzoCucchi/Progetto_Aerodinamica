clear all
close all
tic
%%
%% Input iniziali

% DATI CESNNA 172
b=11;
Delta=0.0000000*pi/180;
d=1*pi/180;

S=16.2;
C_t=1.4;
C_r=S*2/b-C_t;

alfa=0;
beta=0;
[~,~,~,rho]=atmosisa(4100);
U_inf=226/3.6;
c_med=(C_r+C_t)/2;
N=10;
%% Costruisco la geometria
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0],0);
U = U_inf.*[cos(beta*pi/180)*cos(alfa*pi/180) ; -sin(beta*pi/180) ; cos(beta*pi/180)*sin(alfa*pi/180)];

%% Costruisco la geometria dellla coda
% dati cessna 172
S_tail=2;
b_tail=3.4;
C_t_tail=0.4; %DImensione alle estremità
Delta_tail=6/180*pi; % angolo di freccia
d_tail= 0/180*pi ; %Angolo di diedro
C_r_tail=S_tail*2/b_tail-C_t_tail;
c_med_tail=(C_r_tail+C_t_tail)/2;
distance=4+C_r;
hh=-0.2;
% dati 
% dati pamadi
AR=b^2/S;
lamda=C_t/C_r;

lh=distance-C_r/4+C_r_tail/4;
K_AR=1/AR-1/(1+AR^1.7);
K_lamda=(10-3*lamda)/7;
K_r=(1-(hh++distance*sin(4*pi/180))/b)/(2*lh/b)^(1/3);
deps=4.44*(K_AR*K_lamda*K_r*sqrt(cos(Delta)))^1.19
[X_tail,Y_tail,Z_tail,p_tail,f_tail,M_tail]=Geometria(b_tail,Delta_tail,C_r_tail,C_t_tail,d_tail,N,[distance 0 hh],0);
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
%% calcolo a_t e a_w
%%
Q=14;
alfa_vect=linspace(-3,10,Q);
for t=1:Q
    b_tail_s=zeros(f_tail,1);
    b_wing=zeros(f,1);
    b_wing_tail=zeros(f_tot,1);
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
    %calcolo circolazione sulll'ala
    for i=1:f
        Normal=p.panels(i).n;
        b_wing(i)=-dot(U,Normal);
    end
    Gamma_wing=A1\b_wing;
    % Costruisco le matrici
    Gamma_matrix=[];
    for k=0:2*N-1
        Gamma_matrix=[Gamma_matrix;Gamma_alone_tail((1+k*M_tail):(M_tail+k*M_tail))'];
    end
    Gamma_matrix=[fliplr(Gamma_matrix(N+1:end,:)) Gamma_matrix(1:N,:)];
    Gamma_matrix1=[];
    for k=0:2*N-1
        Gamma_matrix1=[Gamma_matrix1;Gamma_wing((1+k*M):(M+k*M))'];
    end
    Gamma_matrix1=[fliplr(Gamma_matrix1(N+1:end,:)) Gamma_matrix1(1:N,:)];
    [Cl_2D_s,L_tail]=portanza(M_tail,Gamma_matrix,rho,U_inf,c_med_tail,d_tail,b_tail,N);
    [Cl_2D_d,L_wing]=portanza(M,Gamma_matrix1,rho,U_inf,c_med,d,b,N);
    Cl_tail_alfa(t)=L_tail/(0.5*rho*U_inf^2*S_tail);
    Cl_wing_alfa(t)=L_wing/(0.5*rho*U_inf^2*S);
    [d_wing_d,gamma_scia]=Resistenza(rho,U_inf,Gamma_matrix1,Y,Z,M,p,alfa,beta,S);
    cd_wing(t)=d_wing_d/(0.5*rho*U_inf^2*S);
    
end

c_t=polyfit(alfa_vect,Cl_tail_alfa,1);
c_w=polyfit(alfa_vect,Cl_wing_alfa,1);
a_w=c_w(1)
a_t=c_t(1)
%% calcolo del modello a due superfici
i_w=4;
i_t=1;
[X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,[0 0 0],i_w);
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
Q=14;
alfa_vect=linspace(-3,10,Q);
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
    [Cl_2D_s,L_tail_d(t)]=portanza(M_tail,Gamma_matrix_d,rho,U_inf,c_med_tail,d_tail,b_tail,N);
    [Cl_2D_d,L_wing_d(t)]=portanza(M,Gamma_matrix1_d,rho,U_inf,c_med,d,b,N);
   % [d_wing_d(t),gamma_scia]=Resistenza(rho,U_inf,Gamma_matrix1_d,Y,Z,M,p,alfa,beta,S);
%[d_tail_d(t),gamma_scia]=Resistenza(rho,U_inf,Gamma_matrix_d,Y_tail,Z_tail,M_tail,p_tail,alfa,beta,S);
  %  Cd(t)=(d_wing_d(t)+d_tail_d(t))/(0.5*rho*U_inf^2*S);
Cl_tot(t)=(L_tail_d(t)+L_wing_d(t))/(0.5*rho*U_inf^2*S);
end
%figure (3)
%plot(Cd,Cl_tot,'b','linewidth',1);

%%
sigma=0.98;
eta=S_tail/S;
for t=1:Q
    eps(t)=i_t-(Cl_tot(t)-i_w*a_w)/(eta*a_t*sigma)+(a_w/eta/a_t/sigma+1)*alfa_vect(t);
    figure (2)
    plot(alfa_vect(t),eps(t),'*','color','r');
    hold on
    
end
c1=polyfit(alfa_vect,eps,1)
for t=1:Q
    eps_p=alfa_vect(t)*c1(1)+c1(2);
CL(t)=a_w*(alfa_vect(t)+i_w)+S_tail/S*a_t*sigma*(alfa_vect(t)+i_t-eps_p);
end
figure (1)
xlrf5=importdata("Progetto aerodinamica1.txt");
xlrf5=xlrf5.data;
plot(alfa_vect,Cl_tot,'b','linewidth',1);
hold on
plot(alfa_vect,CL,'r','linewidth',1);
plot(xlrf5(:,1),xlrf5(:,3),'c','linewidth',1);

grid on
legend("Wessinger","Modello 2 superfici","Xlrf5");
title('Confronto tra i 3 modelli')
xlabel("\alpha °");
ylabel ("CL")
figure (2)
plot(alfa_vect,c1(1).*alfa_vect+c1(2),'b','linewidth',1);
title('Angolo di Downwash in funzione di alfa \epsilon(\alpha):')
xlabel("\alpha °");
ylabel ("\epsilon °")
grid on
CL(6)*(0.5*rho*U_inf^2*S)/9.8;
W=850;
(L_tail_d(6)*(4+C_r_tail/4+C_r)+L_wing_d(6)*(C_r/4))/W;
ee=abs(c1(1)-deps)*2*100/(c1(1)+deps)
%%
err=(xlrf5(:,3)-Cl_tot')./xlrf5(:,3)*100;
%%
CL_3D=load("CL_3D.mat");
CL_3D=CL_3D.CL;
CD_3D=load("CD_3D.mat");
CD_3D=CD_3D.CD;
figure ()
plot(alfa_vect,Cl_wing_alfa,'r','linewidth',1);
hold on
plot(alfa_vect,CL_3D,'b','linewidth',1);
grid on
legend("Weissinger","Linea di Prandtl");
title('Weissinger e Linea di Prandtl: Curva di portanza')
xlabel("\alpha °");
ylabel ("CL")
figure ()
plot(cd_wing,Cl_wing_alfa,'r','linewidth',1);
hold on
plot(CD_3D,CL_3D,'b','linewidth',1);
grid on
legend("Weissinger","Linea di Prandtl");
title('Weissinger e Linea di Prandtl: Curva polare')
xlabel("CD");
ylabel ("CL")
c_prandt=polyfit(alfa_vect,CL_3D,1);
%%
for i=1:14
    err_cd(i)=abs(cd_wing(i)-CD_3D(i))/(cd_wing(i)+CD_3D(i))*2*100;
end
for i=1:14
    err_cl(i)=abs(Cl_wing_alfa(i)-CL_3D(i))/(Cl_wing_alfa(i)+CL_3D(i))*2*100;
end
max(err_cd)
max(err_cl)
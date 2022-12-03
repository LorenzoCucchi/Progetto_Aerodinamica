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
alfa=0;
beta=0;
[~,~,~,rho]=atmosisa(4100);
U_inf=226/3.6;
c_med=(C_r+C_t)/2;
i_w=4;
N=5;
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
distance=4+C_r;
hh=-0.3;
% dati pamadi
AR=b^2/S;
lamda=C_t/C_r;
i_t=1;
c_med_tail=(C_r_tail+C_t_tail)/2
lh=distance-C_r/4+C_r_tail/4;
K_AR=1/AR-1/(1+AR^1.7);
K_lamda=(10-3*lamda)/7;
K_r=(1-hh/b)/(2*lh/b)^(1/3);
deps=4.44*(K_AR*K_lamda*K_r*sqrt(cos(Delta)))^1.19
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
%%
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
%%
Q=11;
alfa_vect=linspace(-5,5,Q);

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
    % calcolo circolazione totale
    for i=1:f_tot
        Normal=p_tot.panels(i).n;
        b_wing_tail(i)=-dot(U,Normal);
    end
    Gamma2=A2\b_wing_tail;
    % separo le due parti
    Gamma_wing_d=Gamma2(1:2*M*N);
    Gamma_tail_d=Gamma2(2*M*N+1:end);
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
    [Cl_2D_s,L_tail]=portanza(M_tail,Gamma_matrix,rho,U_inf,c_med_tail,d_tail,b_tail,N);
    [Cl_2D_d,L_wing]=portanza(M,Gamma_matrix1,rho,U_inf,c_med,d,b,N);
    Cl_tail_alfa(t)=L_tail/(0.5*rho*U_inf^2*S_tail);
    Cl_wing_alfa(t)=L_wing/(0.5*rho*U_inf^2*S);
    % calcolo il Cl totale
    [Cl_2D_s,L_tail_d(t)]=portanza(M_tail,Gamma_matrix_d,rho,U_inf,c_med_tail,d_tail,b_tail,N);
    [Cl_2D_d,L_wing_d(t)]=portanza(M,Gamma_matrix1_d,rho,U_inf,c_med,d,b,N);
    Cl_tot(t)=(L_tail_d(t)+L_wing_d(t))/(0.5*rho*U_inf^2*S);
    
    Cl_S_alfa(t)=L_tail/(0.5*rho*U_inf^2*S_tail);
    Cl_d_alfa(t)=L_tail_d(t)/(0.5*rho*U_inf^2*S_tail);
end
%%
figure ()
plot(alfa_vect,Cl_S_alfa,'b','linewidth',1);
hold on
plot(alfa_vect,Cl_d_alfa,'r','linewidth',1);
grid on
legend("Coda da sola","Coda con ala davanti");
title(['Effetto Downwash:',' i_t',' = ', num2str(i_t),'°',' i_w',' = ', num2str(i_w),'°',' Distanza',' = ', num2str(distance),'m'])
xlabel("\alpha °");
ylabel("Cl");
figure ()

plot(alfa_vect+i_t,Cl_tail_alfa,'b','linewidth',0.5);
hold on
plot(alfa_vect,Cl_wing_alfa,'r','linewidth',0.5);
plot(alfa_vect,Cl_tot,'g','linewidth',0.5);
grid on
legend("Coda","ala ","totale");
title(['Effetto Downwash:',' Calettamento',' = ', num2str(i_t),'°',' Distanza',' = ', num2str(distance),'m'])
xlabel("\alpha");
%%
c1=polyfit(alfa_vect,Cl_tail_alfa,1);
c2=polyfit(alfa_vect,Cl_wing_alfa,1);

m1=c1(1);
m2=c2(1);
i_t=c1(2)/m1
i_w=c2(2)/m2
eta=S_tail/S;
eps=alfa_vect'+i_t.*zeros(Q,1)-(Cl_tot'-m2.*(alfa_vect'+i_w.*zeros(Q,1)))/eta/m1;
c3=polyfit(alfa_vect+i_t.*zeros(1,Q),eps,1);
ee=abs(c3(1)-deps)*2*100/(c3(1)+deps);
fprintf("Errore percentuale: %.2f \n",ee);
%%
G=[C_r/4 0 0 ];
fprintf("Calcolo f_orze.... \n");
b2w=[cos(alfa*pi/180)*cos(beta*pi/180) -sin(beta*pi/180) sin(alfa*pi/180)*cos(beta*pi/180) ;...
    cos(alfa*pi/180)*sin(beta*pi/180)  cos(beta*pi/180) sin(alfa*pi/180)*sin(beta*pi/180) ;...
    -sin(alfa*pi/180)               0         cos(alfa*pi/180)        ];
[F,Moment]=force(p,Gamma,rho,U,G);
F_wind=b2w*F';
CF_wind=F_wind/(0.5*rho*U_inf^2*S);
[F_tail,Moment_tail]=force(p_tail,Gamma_alone_tail,rho,U,G);
F_wind=b2w*F';
%%
CM_wind=[Moment(1)/(0.5*rho*U_inf^2*S*b) Moment(2)/(0.5*rho*U_inf^2*S*c_med) Moment(3)/(0.5*rho*U_inf^2*S*b)];
%Calcolo dei momenti, equibrio ai momenti
%W=1100;
%L_tail_d(6)+L_wing_d(6)-W*9.8

%(L_tail_d(10)*(4+C_r_tail/4+C_r)+L_wing_d(10)*(C_r/4))/W

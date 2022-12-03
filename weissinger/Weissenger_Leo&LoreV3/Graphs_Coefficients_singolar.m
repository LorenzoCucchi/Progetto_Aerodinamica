function Graphs_Coefficients_singolar(p,U_inf,A1,rho,S,b,c_med)
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
    [F_new,M_new]=force(p,Gamma_new,rho,U,C_r);
    F_wind=b2w*F_new';
CF_wind_new(t,:)=F_wind'/(0.5*rho*U_inf^2*S);
CM_wind_new(t,:)=[M_new(1)/(0.5*rho*U_inf^2*S*b) M_new(2)/(0.5*rho*U_inf^2*S*c_med) M_new(3)/(0.5*rho*U_inf^2*S*b)];
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
xlabel("CF");
ylabel("CL");
grid on

%% Calcolo parametri in funzione di Beta 
fprintf("Calcolo polare e curva CL_beta.... \n");
beta_vect=linspace(-5,10,15);
f=2*M*N;
CF_wind_new=zeros(15,3);
CM_wind_new=zeros(15,3);
for t=1:15
    alfa=0;
    beta=beta_vect(t);
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
    [F_new,M_new]=force(p,Gamma_new,rho,U,C_r);
    F_wind=b2w*F_new';
CF_wind_new(t,:)=F_wind'/(0.5*rho*U_inf^2*S);
CM_wind_new(t,:)=[M_new(1)/(0.5*rho*U_inf^2*S*b) M_new(2)/(0.5*rho*U_inf^2*S*c_med) M_new(3)/(0.5*rho*U_inf^2*S*b)];
end
figure()
plot(beta_vect,CF_wind_new(:,3),'r','linewidth',1);
title("Curva CL_\beta");
xlabel("\beta");
ylabel("CL");
grid on
Clm_beta=(CM_wind_new(8,1)-CM_wind_new(7,1))/(beta_vect(8)*pi/180-beta_vect(7)*pi/180);
figure()
plot(CF_wind_new(:,1),CF_wind_new(:,3),'r','linewidth',1);
title("Drag-Polar-beta");
xlabel("CF");
ylabel("CL");
grid on
figure()
plot(beta_vect,CM_wind_new(:,1),'b','linewidth',1);
title("Curva Cl_\beta (Roll)");
xlabel("\beta");
ylabel("Cl_\beta");
grid 
figure()
plot(beta_vect,CM_wind_new(:,3),'b','linewidth',1);
title("Curva Cn_\beta (yaw)");
xlabel("\beta");
ylabel("Cn_\beta");
grid on
end
clear all
close all
clc

%% Richiesta dati del problema

% NACA=input('Inserire il profilo NACA richiesto, a 4 o 5 cifre: ','s');
% N=input('Inserire il numero di pannelli: ');
% alpha=input('Inserire l''angolo d''incidenza in gradi: ');
% U=input('Inserire la velocità di volo in m/s: ');
NACA='0012';
N=160;
alpha=10;
U=1;

alpha=alpha*pi/180;

%% Inizializzazione profilo e pannellizzazione

fprintf('\nInizializzazione profilo...')
[x_p,y_p]=AirfoilShape(NACA,N);
fprintf('\nPannellizzazione...')
[P]=Panels(x_p,y_p);

%% Soluzione del problema

fprintf('\nFormulazione e soluzione del sistema lineare...')
[A]=AICMatrix(P);
[b]=RHS(P,alpha,U);
z=A\b;

fprintf('\nCalcolo delle velocità sul profilo...')
[v]=Velocity(P,alpha,U,z);

fprintf('\nCalcolo del coefficiente di pressione...')
[Cp]=PressureCoeff(v,U);

fprintf('\nCalcolo dei coefficienti di carico...')
[Cl,Cd,Cm]=Loads(P,Cp,U,alpha);

fprintf('\n\nC_l=%0.5g',Cl)

%Calcolo Polare

[Cl_v,Cm_v,Cd_v,alpha_vec] = polar1(P,U);

%% Grafici

fprintf('\n\nGenerazione grafici...\n')

figure(1)
 
subplot(2,1,2) % Profilo
hold on
plot(P.x_v,P.y_v,'k*-') % Vertici dei pannelli
plot(P.x_c,P.y_c,'b*') % Punti medi dei pannelli
quiver(P.x_c,P.y_c,P.n(:,1),P.n(:,2),0.2,'r') % Versori normali
quiver(P.x_v(1:end-1),P.y_v(1:end-1),P.tau(:,1),P.tau(:,2),0.2,'g') % Versori tangenti

axis equal
xlim([-0.1 1.1])
grid on
title("NACA " + NACA + ": profilo",'fontsize',14)
legend('Vertici','Centri','Versori normali','Versori tangenti','interpreter','latex')

subplot(2,1,1) % C_p
hold on
plot(P.x_c(1:N/2),Cp(1:N/2),'r','linewidth',1.5) % Ventre
plot(P.x_c(N/2:end),Cp(N/2:end),'b','linewidth',1.5) % Dorso

set(gca, 'YDir','reverse')
xlim([-0.1 1.1])
grid on
title("NACA " + NACA + ": C_p",'fontsize',14)
legend('Ventre','Dorso','interpreter','latex')


%% Polar Plot

figure(2)
subplot(2,2,1)
hold on
plot(Cd_v,Cl_v,'-^r','LineWidth',1.5)
title("NACA " + NACA + ": CL CD", 'fontsize', 14)

subplot(2,2,2)
plot(alpha_vec*180/pi,Cl_v,'r','linewidth',1.5)
title("NACA " + NACA + ": CL alpha", 'fontsize', 14)

subplot(2,2,3)
plot(alpha_vec*180/pi,Cd_v,'r','LineWidth',1.5)
title("NACA " + NACA + ": CD alpha", 'fontsize', 14)

subplot(2,2,4)
plot(alpha_vec*180/pi, Cm_v,'r','LineWidth',1.5)
title("NACA " + NACA + ": CM alpha", 'fontsize', 14)

















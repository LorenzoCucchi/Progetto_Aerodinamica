clc
clear all
close all

%% general settings

%%% STYLE SETTINGS
set(groot,"defaulttextinterpreter","latex");
set(groot,"defaultAxesTickLabelInterpreter","latex");
set(groot,"defaultLegendInterpreter","latex");
set(groot,"defaultAxesFontSize", 16);
set(0, 'defaultLineLineWidth', 1);
%%% CUSTOM COLORS
colors.red     = [0.6350, 0.0780, 0.1840];
colors.blue    = [0.0000, 0.4470, 0.7410];
colors.orange  = [0.8500, 0.3250, 0.0980];
colors.yellow  = [0.9290, 0.6940, 0.1250];
colors.green   = [0.4660, 0.6740, 0.1880];
colors.azul    = [0.3010, 0.7450, 0.9330];
colors.purple  = [0.4940, 0.1840, 0.5560];
colors.black   = [0.0000, 0.0000, 0.0000];


%plots parameters
N_plot=500;
Rett_calc=true;
Ellisse_calc=false;
Spitfire_calc=true;
Rett2_calc=true;

alfa_inf=-8;
alfa_sup=8;

%ALA RETTANGOLARE CON STESSA CORDA E STESSA SUPERFICIE
if Rett_calc
    %% ALA RETTANGOLARE Cessna 172 SETTINGS
    %general settings
    geom_plots=false;
    distribution_plots=false;
    alfa_plots=false;
    
    
    %field parameters
    Rett.field.Uinf=133.89; % [m/s]
    Rett.field.rho=0.9093; % [kg/m^2]
    
    % Wing parameters
    Rett.wing.b=22.5486/2.54; %apertura alare [m] %dallo spitfire
    
    
    
    Rett.wing.cl_a = 2*pi; %cl derivato alpha del profilo alare (penso debba dipendere da z pure lui ma ok) [1/rad]
    Rett.wing.alfa_g = @(z) deg2rad(0);
    Rett.wing.alfa_0 = @(z) deg2rad(-1.1);
    
    %corda in funzione dell'apertura [m]
    % M=2; %corda massima
    % par_wing.c = @(z) (-M*4/par_wing.b^2)*z.^2 + M ;  %ALA PARABOLICA
    
    Rett.wing.c =@(z)  2.54; %ALA RETTANGOLARE, corda di 1.47 per il Cessna
    
    Rett.wing.S=Rett.wing.b*Rett.wing.c(1); %sup alare [m^2] %16.2 Cessna 172, [modifica il Cl ma non la portanza]
    Rett.wing.AR=Rett.wing.b^2/Rett.wing.S;
    %% Singola configurazione di volo CESSNA 172
    [Rett.Cl, Rett.L_cl, Rett.L_circ, Rett.Gamma, Rett.Cdi, Rett.Di] = single_attitude(Rett.wing,Rett.field,N_plot,distribution_plots);
    
    %% calcolo dei plot rispetto ad alfa (CESSNA 172)
    
    [Rett.Cl_alfa , Rett.Cdi_alfa, Rett.alfa_v] = alpha_plots(Rett.wing,Rett.field,-8,8);
    
    if alfa_plots
        figure()
        plot(Rett.alfa_v,Rett.Cl_alfa)
        title('Cl con alfa')
        grid('minor')
        
        figure()
        hold on
        plot(Rett.alfa_v,Rett.L_cl_alfa)
        
        grid('minor')
        title('Confronto portanze con alfa')
    end

end


%ALA RETTANGOLARE CON STESSA APERTURA E STESSA SUPERFICIE
if Rett2_calc
    %% ALA RETTANGOLARE Cessna 172 SETTINGS
    %general settings
    geom_plots=false;
    distribution_plots=false;
    alfa_plots=false;
    
    
    %field parameters
    Rett2.field.Uinf=133.89; % [m/s]
    Rett2.field.rho=0.9093; % [kg/m^2]
    
    % Wing parameters
    Rett2.wing.b=11.3030; %apertura alare [m]
    
    
    
    Rett2.wing.cl_a = 2*pi; %cl derivato alpha del profilo alare (penso debba dipendere da z pure lui ma ok) [1/rad]
    Rett2.wing.alfa_g = @(z) deg2rad(0);
    Rett2.wing.alfa_0 = @(z) deg2rad(-1.1);
    
    %corda in funzione dell'apertura [m]
    % M=2; %corda massima
    % par_wing.c = @(z) (-M*4/par_wing.b^2)*z.^2 + M ;  %ALA PARABOLICA
    
    Rett2.wing.c =@(z)  22.5486/Rett2.wing.b; %ALA RETTANGOLARE, corda di 1.47 per il Cessna
    
    Rett2.wing.S=Rett2.wing.b*Rett2.wing.c(1); %sup alare [m^2] %16.2 Cessna 172, [modifica il Cl ma non la portanza]
    Rett2.wing.AR=Rett2.wing.b^2/Rett2.wing.S;
    %% Singola configurazione di volo CESSNA 172
    [Rett2.Cl, Rett2.L_cl, Rett2.L_circ, Rett2.Gamma, Rett2.Cdi, Rett2.Di] = single_attitude(Rett2.wing,Rett2.field,N_plot,distribution_plots);
    
    %% calcolo dei plot rispetto ad alfa (CESSNA 172)
    
    [Rett2.Cl_alfa , Rett2.Cdi_alfa, Rett2.alfa_v] = alpha_plots(Rett2.wing,Rett2.field,alfa_inf,alfa_sup);
    
    if alfa_plots
        figure()
        plot(Rett2.alfa_v,Rett2.Cl_alfa)
        title('Cl con alfa')
        grid('minor')
        
        figure()
        hold on
        plot(Rett2.alfa_v,Rett2.L_cl_alfa)

        
        grid('minor')
        title('Confronto portanze con alfa')
    end

end




if Ellisse_calc

    %% ALA ELLITTICA  SETTINGS
    %general settings
    geom_plots=true;
    distribution_plots=false;
    alfa_plots=false;
    
    
    %field parameters
    Ellisse.field.Uinf=51; % [m/s]
    Ellisse.field.rho=1.225; % [kg/m^2]
    
    
    % Wing parameters
    Ellisse.wing.b= 11.3030; %apertura alare [m]
    Ellisse.a=Ellisse.wing.b/2; %semiasse maggiore
    Ellisse.b=1; %semiasse minore
    
   
    Ellisse.wing.S=pi*Ellisse.a*Ellisse.b; %sup alare [m^2] area ellisse
    Ellisse.wing.AR=Ellisse.wing.b^2 / Ellisse.wing.S;
    Ellisse.wing.LAR=1.5; %rapporto longitudinale tra l'ellisse d'uscita e quello d'attacco
    Ellisse.wing.c1=1; %corda composta dall'ellisse d'attacco
    Ellisse.wing.c =@(z) Ellisse.wing.c1*(1 + Ellisse.wing.LAR)*sqrt( Ellisse.b^2 - (z.^2/Ellisse.a^2 * Ellisse.b^2) ); 

    Ellisse.wing.cl_a = 2*pi; %cl derivato alpha del profilo alare (penso debba dipendere da z pure lui ma ok) [1/rad]
    Ellisse.wing.alfa_g = @(z) deg2rad(2);
    Ellisse.wing.alfa_0 = @(z) deg2rad(-1);


   %% Geometry plots
    if geom_plots
        zplot=linspace(-Ellisse.wing.b/2,Ellisse.wing.b/2,500);
        
        figure()
        
        plot(zplot,Ellisse.wing.c(zplot)/(2*Ellisse.wing.c(0)),"Color",'b')
        title('Wing plan')
        hold on
        axis equal
        plot(zplot,-Ellisse.wing.LAR*Ellisse.wing.c(zplot)/(2*Ellisse.wing.c(0)),"Color",'b')
        grid('minor')
      

%         dim=[0.2, 0.25, 0.1, 0.1];
%         annotation('textbox',dim,'String',S_leg,'FitBoxToText','on');
        
        
        
    end
    
    %% Singola configurazione di volo ELLISSE
    [Ellisse.Cl, Ellisse.L_cl, Ellisse.L_circ, Ellisse.Gamma, Ellisse.Cdi, Ellisse.Di] = single_attitude(Ellisse.wing,Ellisse.field,N_plot,distribution_plots);
    
    %% calcolo dei plot rispetto ad alfa ELLISSE
    
    [Ellisse.Cl_alfa , Ellisse.L_cl_alfa , Ellisse.L_circ_alfa , Ellisse.alfa_v] = alpha_plots(Ellisse.wing,Ellisse.field,alfa_inf,alfa_sup);
    
    if alfa_plots
        figure()
        plot(Ellisse.alfa_v,Ellisse.Cl_alfa)
        title('Cl con alfa')
        grid('minor')
        
        figure()
        hold on
        plot(Ellisse.alfa_v,Ellisse.L_cl_alfa)
        plot(Ellisse.alfa_v,Ellisse.L_circ_alfa)
        legend('L_{cl}','L_{circ}')
        
        grid('minor')
        title('Confronto portanze con alfa')
    end

end



if Spitfire_calc

    %% SPITFIRE  SETTINGS
    %general settings
    geom_plots=true;
    distribution_plots=false;
    alfa_plots=false;
    
    
    %field parameters
    Spitfire.field.Uinf=133.89; % [m/s]
    Spitfire.field.rho=0.9093; % [kg/m^2]
    
    
    % Wing parameters

    Spitfire_geom %load Spitfire's data and does the precalculation on the airfoils
    Spitfire.wing.b= 2*span_v(end); %apertura alare [m]
    Spitfire.a=Spitfire.wing.b/2; %semiasse maggiore

    Spitfire.b1=1/4 *chord(1); %semiasse minore dell'ellisse più piccolo
    Spitfire.wing.S1= 0.5 * pi*Spitfire.a*Spitfire.b1; %prima parte della sup alare [m^2] semi-area ellisse più piccolo

    Spitfire.b2=3/4 *chord(1); %semiasse minore dell'ellisse più grande
    Spitfire.wing.S2= 0.5 * pi*Spitfire.a*Spitfire.b2; %seconda parte della sup alare [m^2] semi-area ellisse più grande
    
    Spitfire.wing.S = Spitfire.wing.S1 + Spitfire.wing.S2; %sup alare [m^2] 
    Spitfire.wing.AR=Spitfire.wing.b^2 / Spitfire.wing.S;
    Spitfire.wing.LAR=Spitfire.b2/Spitfire.b1; %rapporto longitudinale tra l'ellisse d'uscita e quello d'attacco
    Spitfire.wing.c1=Spitfire.b1; %corda composta dall'ellisse d'attacco

   
    span_v_interp_tot=[-flip(span_v_interp) , span_v_interp];
    span_v_tot=[-flip(span_v) ; span_v];
    chordfun=polyfit( span_v_tot, [flip(chord) ; chord],14);

    figure()
    hold on
    plot(span_v_tot, [flip(chord) ; chord]);
    axis equal
    plot(span_v_interp_tot,polyval(chordfun,span_v_interp_tot));
    
    Spitfire.wing.c =@(z) polyval(chordfun,z); 

    Spitfire.wing.cl_a = 2*pi; %cl derivato alpha del profilo alare [1/rad]
    Spitfire.wing.alfa_g = @(z) deg2rad(0);
    Spitfire.wing.alfa_0 = @(z) deg2rad(-1.1);


  
    
    
    %% Singola configurazione di volo SPITFIRE
    [Spitfire.Cl, Spitfire.L_cl, Spitfire.L_circ, Spitfire.Gamma, Spitfire.Cdi, Spitfire.Di] = single_attitude(Spitfire.wing,Spitfire.field,N_plot,distribution_plots);
    
    %% calcolo dei plot rispetto ad alfa SPITFIRE
    
    [Spitfire.Cl_alfa , Spitfire.Cdi_alfa , Spitfire.alfa_v] = alpha_plots(Spitfire.wing,Spitfire.field,alfa_inf,alfa_sup);
    
    if alfa_plots
        figure()
        plot(Spitfire.alfa_v,Spitfire.Cl_alfa)
        title('Cl con alfa')

        grid('minor')
        
        figure()
        hold on
        plot(Spitfire.alfa_v,Spitfire.L_cl_alfa)

        
        grid('minor')
        title('Confronto portanze con alfa')
    end

end



%% TOTAL PLOTS

zrettv=linspace(-Rett.wing.b/2,Rett.wing.b/2,N_plot);
zrett2v=linspace(-Rett2.wing.b/2,Rett2.wing.b/2,N_plot);
zspitv=linspace(-Spitfire.wing.b/2,Spitfire.wing.b/2,N_plot);

% DISTRIBUZIONI PORTANZA
figure()
hold on
grid minor

plot(zrettv,-Rett.Gamma*Rett.field.rho*Rett.field.Uinf,'-.',Color=colors.red,LineWidth=1.5)
plot(zrett2v,-Rett2.Gamma*Rett2.field.rho*Rett2.field.Uinf,'--',Color=colors.blue,LineWidth=1.5)
plot(zspitv,-Spitfire.Gamma*Spitfire.field.rho*Spitfire.field.Uinf,Color=colors.green,LineWidth=2)

legend('Rettangolare 1','Rettangolare 2','Spitfire',Location='best')
xlabel('Apertura alare')
ylabel('L')
title('Distribuzione di portanza in apertura alare')
set(gca,'DataAspectRatio',[1 350 1])



% Cl-alfa
% alfa_v=linspace(alfa_inf,alfa_sup,50);
figure()
hold on
grid minor

plot(Rett.alfa_v, Rett.Cl_alfa, '-.',Color=colors.red,LineWidth=1.5)
plot(Rett2.alfa_v, Rett2.Cl_alfa, '--',Color=colors.blue,LineWidth=1.5)
plot(Spitfire.alfa_v, Spitfire.Cl_alfa, Color=colors.green,LineWidth=2)
xline(0)
yline(0)

legend('Rettangolare 1','Rettangolare 2','Spitfire',Location='best')
xlabel('$\alpha$')
ylabel('$C_{l}$')
title('Curva $C_{l}$-$\alpha$')

% Cdi-alfa
figure()
hold on
grid minor

plot(Rett.alfa_v, Rett.Cdi_alfa, '-.',Color=colors.red,LineWidth=1.5)
plot(Rett2.alfa_v, Rett2.Cdi_alfa, '--',Color=colors.blue,LineWidth=1.5)
plot(Spitfire.alfa_v, Spitfire.Cdi_alfa, Color=colors.green,LineWidth=2)
% xline(0)
% yline(0)

legend('Rettangolare 1','Rettangolare 2','Spitfire',Location='best')
xlabel('$\alpha$')
ylabel('$C_{d_{i}}$')
title('Curva $C_{d_{i}}$-$\alpha$')

% Efficienza col Cdi

Rett.E_alfa=abs(Rett.Cl_alfa./Rett.Cdi_alfa);
Rett2.E_alfa=abs(Rett2.Cl_alfa./Rett2.Cdi_alfa);
Spitfire.E_alfa=abs(Spitfire.Cl_alfa./Spitfire.Cdi_alfa);
% 
% figure()
% hold on
% grid minor
% 
% plot(Rett.alfa_v, Rett.E_alfa, '-.',Color=colors.red,LineWidth=1.5)
% plot(Rett2.alfa_v, Rett2.E_alfa, '--',Color=colors.blue,LineWidth=1.5)
% plot(Spitfire.alfa_v, Spitfire.E_alfa, Color=colors.green,LineWidth=2)
% 
% legend('Rettangolare 1','Rettangolare 2','Spitfire',Location='best')
% xlabel('$\alpha$')
% ylabel('$\frac{C_{l}}{C_{d_{i}}}$')
% 
% Rett.E_alfa=Rett.Cl_alfa./Rett.Cdi_alfa;
% Rett2.E_alfa=Rett2.Cl_alfa./Rett2.Cdi_alfa;
% Spitfire.E_alfa=Spitfire.Cl_alfa./Spitfire.Cdi_alfa;
% 

%%%%%%%%%%
figure()
hold on
grid minor
alfa_start=75;
plot(Rett.alfa_v(alfa_start:end), Rett.E_alfa(alfa_start:end), '-.',Color=colors.red,LineWidth=1.5)
plot(Rett2.alfa_v(alfa_start:end), Rett2.E_alfa(alfa_start:end), '--',Color=colors.blue,LineWidth=1.5)
plot(Spitfire.alfa_v(alfa_start:end), Spitfire.E_alfa(alfa_start:end), Color=colors.green,LineWidth=2)

legend('Rettangolare 1','Rettangolare 2','Spitfire',Location='best')
xlabel('$\alpha$')
ylabel('$\frac{C_{l}}{C_{d_{i}}}$')
title('Efficienza aerodinamica calcolata solo con la resistenza indotta')

% Polare col Cdi
figure()
hold on
grid minor

plot(Rett.Cdi_alfa,Rett.Cl_alfa, '-.',Color=colors.red,LineWidth=1.5)
plot(Rett2.Cdi_alfa,Rett2.Cl_alfa, '--',Color=colors.blue,LineWidth=1.5)
plot(Spitfire.Cdi_alfa,Spitfire.Cl_alfa, Color=colors.green,LineWidth=2)

legend('Rettangolare 1','Rettangolare 2','Spitfire',Location='best')
xlabel('$C_{d_{i}}$')
ylabel('$C_{l}$')
title('Curva polare solo con la resistenza indotta')

x=0;
%exportgraphics(gcf,'CN_body.pdf','Resolution',600)
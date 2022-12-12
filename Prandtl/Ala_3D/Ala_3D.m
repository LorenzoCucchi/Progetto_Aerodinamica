clc
clear all
close all

%% general settings

%plots parameters
N_plot=500;
Rett_calc=false;
Ellisse_calc=true;

if Rett_calc
    %% ALA RETTANGOLARE Cessna 172 SETTINGS
    %general settings
    geom_plots=false;
    distribution_plots=true;
    alfa_plots=false;
  
    
    %field parameters
    Rett.field.Uinf=51; % [m/s]
    Rett.field.rho=1.225; % [kg/m^2]
    
    % Wing parameters
    Rett.wing.b=20; %apertura alare [m]
    
    
    
    Rett.wing.cl_a = 6; %cl derivato alpha del profilo alare (penso debba dipendere da z pure lui ma ok) [1/rad]
    Rett.wing.alfa_g = @(z) deg2rad(0);
    Rett.wing.alfa_0 = @(z) deg2rad(-1);
    
    %corda in funzione dell'apertura [m]
    % M=2; %corda massima
    % par_wing.c = @(z) (-M*4/par_wing.b^2)*z.^2 + M ;  %ALA PARABOLICA
    
    Rett.wing.c =@(z)  2.6667; %ALA RETTANGOLARE, corda di 1.47 per il Cessna
    
    Rett.wing.S=Rett.wing.b*Rett.wing.c(1); %sup alare [m^2] %16.2 Cessna 172, [modifica il Cl ma non la portanza]
    Rett.wing.AR=Rett.wing.b^2/Rett.wing.S;
    %% Singola configurazione di volo CESSNA 172
    [Rett.Cl, Rett.L_cl, Rett.L_circ, Rett.Gamma, Rett.Cdi] = single_attitude(Rett.wing,Rett.field,N_plot,distribution_plots);
    
    %% calcolo dei plot rispetto ad alfa (CESSNA 172)
    
    [Rett.Cl_alfa , Rett.L_cl_alfa , Rett.L_circ_alfa , Rett.alfa_v] = alpha_plots(Rett.wing,Rett.field,-5,12);
    
    if alfa_plots
        figure()
        plot(Rett.alfa_v,Rett.Cl_alfa)
        title('Cl con alfa')
        grid('minor')
        
        figure()
        hold on
        plot(Rett.alfa_v,Rett.L_cl_alfa)
        plot(Rett.alfa_v,Rett.L_circ_alfa)
        legend('L_{cl}','L_{circ}')
        
        grid('minor')
        title('Confronto portanze con alfa')
    end

end

if Ellisse_calc

    %% ALA ELLITTICA  SETTINGS
    %general settings
    geom_plots=true;
    distribution_plots=true;
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

    Ellisse.wing.cl_a = 6; %cl derivato alpha del profilo alare (penso debba dipendere da z pure lui ma ok) [1/rad]
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
    
    [Ellisse.Cl_alfa , Ellisse.L_cl_alfa , Ellisse.L_circ_alfa , Ellisse.alfa_v] = alpha_plots(Ellisse.wing,Ellisse.field,-5,12);
    
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
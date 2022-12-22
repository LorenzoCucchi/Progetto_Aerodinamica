%% NACA FINDER
clc
clear all
close all


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



super_cl_alfa=cl_alfa_import("supermarine371I.dat");
super_polar=polar_import("supermarine371I.dat");

super.cl_alfa=super_cl_alfa;
super.polar=super_polar;
% 
% figure()
% plot(naca(:,1),naca(:,2))

% af='2313';
% alpha_and_t=13:-1:6;

alpha=linspace(-15,15,100);
Re=1e6;
Mach=0;

naca_des=[];
t="13"; %1413 1206 [1413 1412 1411 1410 1209 1208 1207 1206]
MP=1:1:20;
% MP=[MP, 45:1:49]; %, 55:1:59, 65:1:69, 75:1:79, 85:1:89, 95:1:99
for n=1:length(MP)
    mp=MP(n);
    if(mp<10)
        nome=strcat("0",string(mp));
        naca_des=[strcat(nome,t); naca_des];
        
    else
        naca_des=[strcat(string(mp),t); naca_des];
    end
end


[af_found, diff_found] = find_airfoil(naca_des,alpha,super); %ROOT: 1413, TIP: 1007
% af_found= "1207";
Re=1e6;
Mach=0;
[pol_found,foil_found] = xfoil(convertStringsToChars(strcat("NACA ",af_found)),alpha,Re,Mach);



%CL-ALFA
figure()
hold on 
plot(super_cl_alfa(:,1),super_cl_alfa(:,2),'Color',[0 0.4470 0.7410])
plot(pol_found.alpha,pol_found.CL,'Color','black')
legend('Supermarine','Naca')
xlabel('$\alpha$')
ylabel('$C_{l}$')
grid minor

%DIFF CL-ALFA
figure()
hold on 
plot(alpha,diff_found)
title('Difference in CL from the objective')
grid minor

%POLAR
figure()
hold on 
plot(super_polar(:,2),super_polar(:,1),'Color',[0 0.4470 0.7410])
plot(pol_found.CD,pol_found.CL,'Color','black')
grid minor
legend('Supermarine','Naca')

%GEOMETRY
airfoil.designation=convertStringsToChars(af_found);
airfoil.n=50;
airfoil.HalfCosineSpacing=1;
airfoil.wantFile=0;
airfoil.is_finiteTE=0;
plot_af_found = naca4gen(airfoil);

figure()
hold on
load("afsupermarineI.mat")
plot(afsupermarineI(:,1),afsupermarineI(:,2),'Color',[0 0.4470 0.7410])
plot(plot_af_found.x,plot_af_found.z,'Color','black')
legend('Supermarine','Naca')
axis equal

%exportgraphics(gcf,'CN_body.pdf','Resolution',600)
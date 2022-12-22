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



single_plot_flag=true;
alpha=linspace(-3,3,150);
Re=1e8;
Mach=0;

%a0=[-1.0796	-1.0182	-1.0007 -0.9941 -0.9139 -0.9146 -0.9150 -0.8794]
naca_des=[];
t_v=["06" "07" "08" "09" "10" "11" "12" "13"]; %1413 1206 [1413 1412 1411 1410 1209 1208 1207 1206]
% t_v=["07" "06"];
MP=14;
ao=nan(length(t_v),1);
for k=1:length(t_v)
    
%     if str2double(t_v(k))>8
%         alpha=linspace(-9,9,120);
%     end
    if(MP<10)
        nome=strcat("0",string(MP));
        naca_des=strcat(nome,t_v(k));
        
    else
        naca_des=strcat(string(MP),t_v(k));
    end
     
    [pol,foil] = xfoil(convertStringsToChars(strcat("NACA ",naca_des)),alpha,Re,Mach);
    
    cl_interp=polyfit(pol.alpha,pol.CL,1);
    cl_a=rad2deg(polyder(cl_interp));
    a0(k)=roots(cl_interp);

    
   
    if single_plot_flag
        
        %CL-ALFA
%         figure()
%         hold on 
%         plot(pol.alpha,pol.CL,'Color','black')
%         grid minor
%         
%         %POLAR
%         figure()
%         hold on 
%         plot(pol.CD,pol.CL,'Color','black')
%         grid minor
%         
%         
%         %GEOMETRY
%         airfoil.designation=convertStringsToChars(naca_des);
%         airfoil.n=50;
%         airfoil.HalfCosineSpacing=1;
%         airfoil.wantFile=0;
%         airfoil.is_finiteTE=0;
%         plot_af = naca4gen(airfoil);
%         
%         figure()
%         hold on
%         plot(plot_af.x,plot_af.z,'Color','black')
%         axis equal
        
        %CL-alfa slope calculation
        
        figure()
        hold on
        plot(pol.alpha,pol.CL,'black')
        plot(pol.alpha,polyval(cl_interp,pol.alpha),'red')
        grid minor
        title('interpolated CL-alpha curve')
    end
end

figure()
plot(str2double(t_v),a0)
grid minor
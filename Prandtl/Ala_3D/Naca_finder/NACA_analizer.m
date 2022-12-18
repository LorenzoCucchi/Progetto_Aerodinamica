clc
clear all
close all

plot_flag=false;
alpha=linspace(-5,5,30);
Re=1e6;
Mach=0;

%a0=[-1.0796	-1.0182	-1.0007 -0.9941 -0.9139 -0.9146 -0.9150 -0.8794]
naca_des=[];
t_v=["07" "06"]; %1413 1206 [1413 1412 1411 1410 1209 1208 1207 1206]
MP=12;
ao=nan(length(t_v),1);
for k=1:length(t_v)
    
        
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

    if plot_flag
        
        %CL-ALFA
        figure()
        hold on 
        plot(pol.alpha,pol.CL,'Color','black')
        grid minor
        
        %POLAR
        figure()
        hold on 
        plot(pol.CD,pol.CL,'Color','black')
        grid minor
        
        
        %GEOMETRY
        airfoil.designation=convertStringsToChars(naca_des);
        airfoil.n=50;
        airfoil.HalfCosineSpacing=1;
        airfoil.wantFile=0;
        airfoil.is_finiteTE=0;
        plot_af = naca4gen(airfoil);
        
        figure()
        hold on
        plot(plot_af.x,plot_af.z,'Color','black')
        axis equal
        
        %CL-alfa slope calculation
        
        figure()
        hold on
        plot(pol.alpha,pol.CL,'black')
        plot(pol.alpha,polyval(cl_interp,pol.alpha),'red')
        grid minor
        title('interpolated CL-alpha curve')
    end
end
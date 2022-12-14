% thickness distribution along the wingspan of the Supermarine Spitfire
%NACA MPXX
%M=max camber % 
%P=position of the maximum camber divided by 10
%XX=max thickness % (T)

%% import cl_alfa
clc
clear all
close all


super_cl_alfa=cl_alfa_import("supermarine371.dat");
super_polar=polar_import("supermarine371.dat");

super.cl_alfa=super_cl_alfa;
super.polar=super_polar;
% 
% figure()
% plot(naca(:,1),naca(:,2))

af='2313';
alpha=linspace(-12,12,30);
Re=1e6;
Mach=0;

naca_des=[];
t="13";
for n=1:40
    if(n<10)
        nome=strcat("0",string(n));
        naca_des=[strcat(nome,t); naca_des];
        
    else
        naca_des=[strcat(string(n),t); naca_des];
    end
end


af_found = find_airfoil(naca_des,alpha,super);
Re=1e6;
Mach=0;
[pol_found,foil_found] = xfoil(convertStringsToChars(strcat("NACA ",af_found)),alpha,Re,Mach);

%CL-ALFA
figure()
hold on 
plot(super_cl_alfa(:,1),super_cl_alfa(:,2))
plot(pol_found.alpha,pol_found.CL,'Color','black')
legend('Supermarine','Naca')

%POLAR
figure()
hold on 
plot(super_polar(:,2),super_polar(:,1))
plot(pol_found.CD,pol_found.CL,'Color','black')

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
load("afsupermarine.mat")
plot(afsupermarine(:,1),afsupermarine(:,2))
plot(plot_af_found.x,plot_af_found.z,'Color','black')
axis equal
%% spitfire
load('span_thickness.mat');
load('chord.mat');

span_v=span_thickness(:,1);
T_v=span_thickness(:,2)./chord *100; %thickness % over chord length (XX)

%from tables:
T_v(end-2)=7.25;
T_v(end-1)=6.72;
T_v(end)=6.15;
figure()
plot(span_v,T_v)
grid minor
ylabel('maximum airfoil thickness %')
xlabel('wingspan')
title('maximum thickness of airfoil through the wingspan')


% T_v contains XX
% set P=40
%set M= 2
M='2';
P='8';
formatSpec = '%.1f';
T_str1=num2str(T_v(T_v>10),formatSpec);
T_str2=strcat('0',num2str(T_v(T_v<10),formatSpec))
T_str=[T_str1;
        T_str2]
NACAs=strcat(M,P,T_str)



for i=1:length(NACAs)

    airfoil.designation=NACAs(i,1:6);
    airfoil.n=50;
    airfoil.HalfCosineSpacing=1;
    airfoil.wantFile=0;
    airfoil.is_finiteTE=0;
    af(i) = naca4gen(airfoil);

end


figure()
hold on
for i=1:24
    plot(af(i).x,af(i).z,'Color','black')
end
axis equal
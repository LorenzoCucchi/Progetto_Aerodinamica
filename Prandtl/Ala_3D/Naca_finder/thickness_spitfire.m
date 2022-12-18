% thickness distribution along the wingspan of the Supermarine Spitfire
%NACA MPXX
%M=max camber % 
%P=position of the maximum camber divided by 10
%XX=max thickness % (T)

%% NACA FINDER
clc
clear all
close all


super_cl_alfa=cl_alfa_import("supermarine371II.dat");
super_polar=polar_import("supermarine371II.dat");

super.cl_alfa=super_cl_alfa;
super.polar=super_polar;
% 
% figure()
% plot(naca(:,1),naca(:,2))

% af='2313';
% alpha_and_t=13:-1:6;

alpha=linspace(-3.5,3.5,30);
Re=1e6;
Mach=0;

naca_des=[];
t="06"; %1413 1206 [1413 1412 1411 1410 1209 1208 1207 1206]
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
load("afsupermarineII.mat")
plot(afsupermarineII(:,1),afsupermarineII(:,2),'Color',[0 0.4470 0.7410])
plot(plot_af_found.x,plot_af_found.z,'Color','black')
legend('Supermarine','Naca')
axis equal

%CL-alfa slope calculation
cl_interp=polyfit(pol_found.alpha,pol_found.CL,1);
cl_a=rad2deg(polyder(cl_interp));
a0=
figure()
hold on
plot(pol_found.alpha,pol_found.CL,'black')
plot(pol_found.alpha,polyval(cl_interp,pol_found.alpha),'red')
grid minor
title('interpolated CL-alpha curve')
%% spitfire
clear all
close all
clc

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
M='1';
P='4';
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

%NACA airfoils
figure()
hold on
for i=1:24
    plot(af(i).x,af(i).z,'Color','black')
end
axis equal
title('Wing discretized into 24 airfoils')

%Airfoils comparison
figure()
hold on

%real
load("afsupermarineI.mat")
plot(afsupermarineI(:,1),afsupermarineI(:,2),'--','Color',[0 0.4470 0.7410],'LineWidth',1.5)
load("afsupermarineII.mat")
plot(afsupermarineII(:,1),afsupermarineII(:,2),'--','Color',[0 0.4470 0.7410])

%naca
plot(af(1).x,af(1).z,'Color','black','LineWidth',1.5)
plot(af(end).x,af(end).z,'Color','black')

axis equal
title('Root and tip airfoils')
legend('Root Supermarine','Tip Supermarine','Root NACA','Tip NACA')


%3D wing
figure()
hold on
axis equal
view(3)

x=nan(length(af(1).x),24);
y=nan(length(af(1).x),24);
z=nan(length(af(1).x),24);
for i=1:24

    x(:,i)=-chord(i)/4 + chord(i)*af(i).x;
    y(:,i)=span_v(i)*ones(length(af(i).x),1);
    z(:,i)=chord(i)*af(i).z;
    plot3(x,y,z,'Color','none')
    
    if i>1
        s=surf([x(:,i-1)'; x(:,i)'], [y(:,i-1)'; y(:,i)'], [z(:,i-1)'; z(:,i)'], 'MeshStyle','row');    
%         s.FaceColor="#77AC30";
        s.FaceColor='interp';
        s.EdgeColor="#77AC30";
    end
end



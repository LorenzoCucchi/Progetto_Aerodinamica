% thickness distribution along the wingspan of the Supermarine Spitfire
%NACA MPXX
%M=max camber % 
%P=position of the maximum camber divided by 10
%XX=max thickness % (T)
clc
clear all
close all

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


figure()
hold on
for i=1:24
    plot(af(i).x,af(i).z,'Color','black')
end
axis equal
% thickness distribution along the wingspan of the Supermarine Spitfire
%NACA MPXX
%M=max camber % 
%P=position of the maximum camber divided by 10
%XX=max thickness % (T)

%% spitfire
% clear all
% close all
% clc

load('span_thickness.mat');
load('chord.mat');
load('a0_14.mat'); %contiene gli a0 dei profili da 1406 a 1413

chord(end)=0.3;
span_v=span_thickness(:,1);
T_v=span_thickness(:,2)./chord *100; %thickness % over chord length (XX)

%from tables:
T_v(end-2)=7.25;
T_v(end-1)=6.72;
T_v(end)=6.15;
figure()
plot(span_v,T_v)
grid minor
ylabel('Massimo spessore del profilo alare \%')
xlabel('Apertura alare')
title('Massimo spessore del profilo lungo l''apertura alare')


% T_v contains XX
% set P=40
%set M= 2
M='1';
P='1';
formatSpec = '%.1f';
T_str1=num2str(T_v(T_v>10),formatSpec);
T_str2=strcat('0',num2str(T_v(T_v<10),formatSpec));
T_str=[T_str1;
        T_str2];
NACAs=strcat(M,P,T_str);



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
title('Normalized airfoils used for discretization')

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
title('Profili alari alla radice e all''estremit\`a')
legend('Radice Supermarine','Estremit\`a Supermarine','Radice NACA','Estremit\`a NACA')



%3D wing
figure()
hold on
grid on
axis equal
xlim([-2 2.5])
ylim([-10 10])
zlim([-1.5 1.5])
view(3)
axis vis3d

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
        s.FaceColor="#77AC30";
%         s.FaceColor='interp';
%         s.EdgeColor="#77AC30";
        s.EdgeColor="black";
    end
end

x=nan(length(af(1).x),24);
y=nan(length(af(1).x),24);
z=nan(length(af(1).x),24);
for i=1:24

    x(:,i)=-chord(i)/4 + chord(i)*af(i).x;
    y(:,i)=-span_v(i)*ones(length(af(i).x),1);
    z(:,i)=chord(i)*af(i).z;
    plot3(x,y,z,'Color','none')
    
    if i>1
        s=surf([x(:,i-1)'; x(:,i)'], [y(:,i-1)'; y(:,i)'], [z(:,i-1)'; z(:,i)'], 'MeshStyle','row');    
        s.FaceColor="#77AC30";
%         s.FaceColor='interp';
        s.EdgeColor="black";
    end
end




%a0 with wing span
%si può fare automaticamente ma non so come
T_str_string=["6" "6" "7" "7" "8" "8" "9" "9" "9" "10" "10" "10" "10" "11" "11" "11" "11" "11" "12" "12" "12" "12" "12" "13"];

a0_span_naca=nan(length(T_str_string),1);

for i=1:length(T_str_string)

    if(str2double(T_str_string(i))==6)
        a0_span_naca(i)=a0(1);
    elseif(str2double(T_str_string(i))==7)
        a0_span_naca(i)=a0(2);

    elseif(str2double(T_str_string(i))==8)
        a0_span_naca(i)=a0(3);

    elseif(str2double(T_str_string(i))==9)
        a0_span_naca(i)=a0(4);

    elseif(str2double(T_str_string(i))==10)
        a0_span_naca(i)=a0(5);

    elseif(str2double(T_str_string(i))==11)
        a0_span_naca(i)=a0(6);

    elseif(str2double(T_str_string(i))==12)
        a0_span_naca(i)=a0(7);

    elseif(str2double(T_str_string(i))==13)
        a0_span_naca(i)=a0(8);
    end

end

a0_span_interp=polyfit(span_v,a0_span_naca,4);
span_v_interp=linspace(span_v(1),span_v(end),100);


figure()
hold on
plot(span_v,a0_span_naca)
plot(span_v_interp,polyval(a0_span_interp,span_v_interp))
grid minor
title('alpha 0 with the wingspan')
legend('discretized','fitted')
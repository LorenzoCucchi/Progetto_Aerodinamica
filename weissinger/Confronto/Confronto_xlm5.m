clear all
close all
Fr_30=load("Fr_30.mat");
Fr_20=load("Fr_20.mat");
Fr_10=load("Fr_10.mat");
Fr_0=load("Fr_0.mat");
d_0=load("d_0.mat");
d_8=load("d_8.mat");
b_10=load("b_10.mat");
b_8=load("b_8.mat");
b_6=load("b_6.mat");
Standard_xlf5 = importdata('Ala standard (d4,b8,Fr0).txt');
% Fr_30_xlf5.data.alfa(1) Cl(3) Cd(4) cm(9)
Fr_30_xlf5=importdata("Fr30.txt");
Fr_20_xlf5=importdata("Fr20.txt");
Fr_10_xlf5=importdata("Fr10.txt");
d_0_xlf5=importdata("d0.txt");
d_8_xlf5=importdata("d8.txt");
b_10_xlf5=importdata("b10.txt");
b_8_xlf5=importdata("b8.txt");
b_6_xlf5=importdata("b6.txt");
Standard_xlf5=Standard_xlf5.data;
Fr_30_xlf5=Fr_30_xlf5.data;
Fr_20_xlf5=Fr_20_xlf5.data;
Fr_10_xlf5=Fr_10_xlf5.data;
d_0_xlf5=d_0_xlf5.data;
d_8_xlf5=d_8_xlf5.data;
b_10_xlf5=b_10_xlf5.data;
b_8_xlf5=b_8_xlf5.data;
b_6_xlf5=b_6_xlf5.data;
%% Variazione angolo di freccia
%% calcolo degli errori in percentuale
err_CL=[];
for i=1:10
    b(i)=Fr_0.ris.num(1+i).CL;
end
err_CL_i=abs(b'-Standard_xlf5(:,3))./(b'+Standard_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
for i=1:10
    b(i)=Fr_10.ris.num(1+i).CL;
end
err_CL_i=abs(b'-Fr_10_xlf5(:,3))./(b'+Fr_10_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
for i=1:10
    b(i)=Fr_20.ris.num(1+i).CL;
end
err_CL_i=abs(b'-Fr_20_xlf5(:,3))./(b'+Fr_20_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
for i=1:10
    b(i)=Fr_30.ris.num(1+i).CL;
end
err_CL_i=abs(b'-Fr_30_xlf5(:,3))./(b'+Fr_30_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
for i=1:11
    b(i)=b_8.ris.num(i).CL;
end
err_CL_i=abs(b'-b_8_xlf5(:,3))./(b'+b_8_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
clear b
for i=1:10
    b(i)=b_10.ris.num(1+i).CL;
end
err_CL_i=abs(b'-b_10_xlf5(:,3))./(b'+b_10_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
for i=1:10
    b(i)=b_6.ris.num(1+i).CL;
end
err_CL_i=abs(b'-b_6_xlf5(:,3))./(b'+b_6_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
for i=1:10
    b(i)=d_8.ris.num(1+i).CL;
end
err_CL_i=abs(b'-d_8_xlf5(:,3))./(b'+d_8_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
for i=1:10
    b(i)=d_0.ris.num(1+i).CL;
end
err_CL_i=abs(b'-d_0_xlf5(:,3))./(b'+d_0_xlf5(:,3)).*2;
err_CL=[err_CL;err_CL_i];
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
h=histogram(err_CL*100,5,'BinWidth',0.5,'Facealpha',0.9,'Facecolor','r')
h.NumBins=10;
ylabel("Count")
xlabel("Err%")
grid on
title("Differenza percentuale tra il codice e xflr5 per CL ")
%% errore cD
err_CD=[];
for i=1:10
    b(i)=Fr_0.ris.num(1+i).CD;
end
err_CD_i=abs(b'-Standard_xlf5(:,4))./(b'+Standard_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
size(err_CD,1)

for i=1:10
    b(i)=Fr_10.ris.num(1+i).CD;
end
err_CD_i=abs(b'-Fr_10_xlf5(:,4))./(b'+Fr_10_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
size(err_CD,1)
for i=1:10
    b(i)=Fr_20.ris.num(1+i).CD;
end
err_CD_i=abs(b'-Fr_20_xlf5(:,4))./(b'+Fr_20_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
size(err_CD,1)

for i=1:10
    b(i)=Fr_30.ris.num(1+i).CD;
end
err_CD_i=abs(b'-Fr_30_xlf5(:,4))./(b'+Fr_30_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
for i=1:11
    b(i)=b_8.ris.num(i).CD;
end

err_CD_i=abs(b'-b_8_xlf5(:,4))./(b'+b_8_xlf5(:,4)).*2;
clear b
err_CD=[err_CD;err_CD_i];
for i=1:10
    b(i)=b_10.ris.num(1+i).CD;
end
err_CD_i=abs(b'-b_10_xlf5(:,4))./(b'+b_10_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
for i=1:10
    b(i)=b_6.ris.num(1+i).CD;
end
err_CD_i=abs(b'-b_6_xlf5(:,4))./(b'+b_6_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
for i=1:10
    b(i)=d_8.ris.num(1+i).CD;
end
err_CD_i=abs(b'-d_8_xlf5(:,4))./(b'+d_8_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
for i=1:10
    b(i)=d_0.ris.num(1+i).CD;
end
err_CD_i=abs(b'-d_0_xlf5(:,4))./(b'+d_0_xlf5(:,4)).*2;
err_CD=[err_CD;err_CD_i];
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
h=histogram(err_CD*100,10,'BinWidth',0.5,'Facealpha',0.9,'Facecolor','r')
h.NumBins=10;
ylabel("Count")
xlabel("Err%")
grid on
title("Differenza percentuale tra il codice e xflr5 per CD ")
%% errore cM
err_CM=[];
for i=1:10
    b(i)=Fr_0.ris.num(1+i).CMa;
end
err_CMi=abs(b'-Standard_xlf5(:,9))./(b'+Standard_xlf5(:,9)).*2;
err_CM=[err_CM;err_CMi];
size(err_CM,1)

for i=1:10
    b(i)=Fr_10.ris.num(1+i).CMa;
end
err_CMi=abs(b'-Fr_10_xlf5(:,9))./(b'+Fr_10_xlf5(:,9)).*2;
err_CM=[err_CM;err_CMi];
size(err_CM,1)
for i=1:10
    b(i)=Fr_20.ris.num(1+i).CMa;
end
err_CMi=abs(b'-Fr_20_xlf5(:,9))./(b'+Fr_20_xlf5(:,9)).*2;
err_CM=[err_CM;err_CMi];
size(err_CM,1)

for i=1:10
    b(i)=Fr_30.ris.num(1+i).CMa;
end
err_CMi=abs(b'-Fr_30_xlf5(:,9))./(b'+Fr_30_xlf5(:,9)).*2;
err_CM=[err_CM;err_CMi];
for i=1:11
    b(i)=b_8.ris.num(i).CMa;
end

err_CMi=abs(b'-b_8_xlf5(:,9))./(b'+b_8_xlf5(:,9)).*2;
clear b
err_CM=[err_CM;err_CMi];
for i=1:10
    b(i)=b_10.ris.num(1+i).CMa;
end
err_CMi=abs(b'-b_10_xlf5(:,9))./(b'+b_10_xlf5(:,9)).*2;
err_CM=[err_CM;err_CMi];
for i=1:10
    b(i)=b_6.ris.num(1+i).CMa;
end
err_CMi=abs(b'-b_6_xlf5(:,9))./(b'+b_6_xlf5(:,9)).*2;
err_CM=[err_CM;err_CMi];
for i=1:10
    b(i)=d_8.ris.num(1+i).CMa;
end
err_CMi=abs(b'-d_8_xlf5(:,9))./(b'+d_8_xlf5(:,9)).*2;
err_CM=[err_CM;err_CMi];
for i=1:10
    b(i)=d_0.ris.num(1+i).CMa;
end
err_CMi=abs(b'-d_0_xlf5(:,9))./(b'+d_0_xlf5(:,9)).*2;
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
err_CM=[-err_CM;-err_CMi];
h=histogram(err_CM*100,10,'BinWidth',0.5,'Facealpha',0.9,'Facecolor','r')
h.NumBins=10;
ylabel("Count")
xlabel("Err%")
grid on
title("Differenza percentuale tra il codice e xflr5 per CM_a ")
%%
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
diedro=30;
d=4;
span=7.905;
for i=1:10
    a(i)=Fr_30.ris.num(1+i).alfa;
    b(i)=Fr_30.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','b');
hold on
fig=plot(Fr_30_xlf5(:,1),Fr_30_xlf5(:,3),'linewidth',0.5,'color','r');

legend("Codice","xflr5");
title(['Curva CL_\alpha Freccia:',' = ', num2str(diedro),'°',' diedro:',' = ', num2str(d),'°',' b',' = ', num2str(span),'m'])
xlabel("\alpha");
ylabel("CL");
grid on
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
diedro=30;
d=4;
span=7.905;
for i=1:10
    a(i)=Fr_30.ris.num(1+i).CD;
    b(i)=Fr_30.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','b');
hold on
plot(Fr_30_xlf5(:,4),Fr_30_xlf5(:,3),'linewidth',0.5,'color','r');
legend("Codice","xflr5");
title(['Curva polare Freccia:',' = ', num2str(diedro),'°',' diedro:',' = ', num2str(d),'°',' b',' = ', num2str(span),'m'])
xlabel("CD");
ylabel("CL");
grid on
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
diedro=30;
d=4;
span=7.905;
for i=1:10
    a(i)=Fr_30.ris.num(1+i).alfa;
    b(i)=Fr_30.ris.num(1+i).CMa;
end
plot(a,b,'linewidth',1,'color','b');
hold on
plot(Fr_30_xlf5(:,1),Fr_30_xlf5(:,9),'linewidth',0.5,'color','r');
legend("Codice","xflr5");
title(['Curva CM_\alpha Freccia:',' = ', num2str(diedro),'°',' diedro:',' = ', num2str(d),'°',' b',' = ', num2str(span),'m'])
xlabel("\alpha");
ylabel("CM");
grid on
%% Variazione
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
d=4;
span=7.905;
for i=1:10
    a(i)=Fr_0.ris.num(1+i).CD;
    b(i)=Fr_0.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','r');
hold on
for i=1:10
    a(i)=Fr_10.ris.num(1+i).CD;
    b(i)=Fr_10.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','m');
for i=1:10
    a(i)=Fr_20.ris.num(1+i).CD;
    b(i)=Fr_20.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','c');
for i=1:10
    a(i)=Fr_30.ris.num(1+i).CD;
    b(i)=Fr_30.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','g');
legend("\Lambda=0°","\Lambda=10°","\Lambda=20°","\Lambda=30°");
title("Curva polare al variare dell'angolo di Freccia \Lambda")
xlabel("CL");
ylabel("CD");
grid on
%% Variazione
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
d=4;
span=7.905;
for i=1:10
    a(i)=b_6.ris.num(1+i).CD;
    b(i)=b_6.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','r');
hold on
for i=1:10
    a(i)=b_8.ris.num(1+i).CD;
    b(i)=b_8.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','m');
for i=1:10
    a(i)=b_10.ris.num(1+i).CD;
    b(i)=b_10.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','g');
legend("b=6m","b=8m","b=10m");
title("Curva CL_\alpha al variare dell'apertura alare b")
xlabel("\alpha");
ylabel("CL");
grid on
%% Variazione
figure ()
set(gcf, 'Position',  [100, 100, 500, 400])
d=4;
span=7.905;
for i=1:10
    a(i)=d_0.ris.num(1+i).CD;
    b(i)=d_0.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','r');
hold on
for i=1:10
    a(i)=Fr_0.ris.num(1+i).CD;
    b(i)=Fr_0.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','r');
hold on
for i=1:10
    a(i)=d_8.ris.num(1+i).CD;
    b(i)=d_8.ris.num(1+i).CL;
end
plot(a,b,'linewidth',1,'color','m');

legend("d=0°","d=4°","d=8°","d=12°");
title("Curva CL_\alpha al variare dell'angolo di diedro b")
xlabel("\alpha");
ylabel("CL");
grid on
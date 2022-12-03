function [X,Y,Z,p,f,M]=Geometria(b,Delta,C_r,C_t,d,N,P_0,theta)
% La function permette di generare tutti i dati necessari per inizializzare
% la funzione
% INPUT
% - b:      apertura alare
% - Delta:  angolo di freccia
% - C_r:    lunghezza della corda alla radice
% - C_t:    lunghezza della corda all'estremità alare
% - d:      angolo di diero
% - N:      numero di suddivisioni lungo la corda
% OUTPUT
% - X:      Matrice (n+1)X(m+1): X dei punti estremi del pannello
% - Y:      Matrice (n+1)X(m+1): Y dei punti estremi del pannello
% - Z:      Matrice (n+1)X(m+1): Z dei punti estremi del pannello
% - P:      Struct: Contiene tutte le informazioni del pannello
%           -P1,P2,P3,P4: estremi del pannello [1x3]
%           -A,B,C,D: punti del vortice a forma di ferro di cavallo
%           -P: Punto di controllo [1x3]
%           -n: normale al pannello [1x3]
% - f:      Numero di pannelli
% - M:      fattore di pannellizzazione in direzione y (semiala)
% Definisco alcuni parametri
c_med=(C_r+C_t)/2;
P_0c=P_0+[C_r/4 0 0];
%M=round(N*b/2/c_med,0);
M=15;
theta=-theta*pi/180;
T=[cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)];
% Adimensionalizzo
%C_r=C_r/c_med;
%C_t=C_t/c_med;
%b=b/c_med;
% Definiscono le funzioni per il calcolo delle coordinate
fprintf("Ricerca dei Punti di Estremità.... \n");
X_r=@(i) C_r./N.*i;
X_t=@(i) C_t./N.*i;
tan_alfa=@(i) (b./2.*tan(Delta)+X_t(i)-X_r(i)).*2./b;
X_pij=@(i,j) P_0(1)+X_r(i)+b./2./M.*j.*tan_alfa(i);
Y_pij=@(j) P_0(2)+b./2./M.*j;
Z_pij=@(j) P_0(3)+b./2./M.*j.*tan(d);
A=zeros((N+1)*(M+1),3);
X=[];
Y=[];
Z=[];
f=1;
for i=0:N
    for j=0:M
        A(f,:)=(T*([X_pij(i,j) Y_pij(j) Z_pij(j)]-P_0c)')'+P_0c;
        X(i+1,j+1)= A(f,1);
        Y(i+1,j+1)= A(f,2);
        Z(i+1,j+1)= A(f,3);
        f=f+1;
    end
end
A=[A;A(:,1) -A(:,2) A(:,3)];
%Calcolo la matrice p per calcolare tutti i pannelli con le varie
%informazioni
p=struct;
f=1;
t=1;
i=1;
fprintf("Creazione dei pannelli...\n");
while i<=((N+1)*(M+1)-M-2)
    %Controllo se sono arrivato al bordo
    if i/(M+1)==t
        t=t+1;
        i=i+1;
    end
    p.panels(f).P4=P_0c+(T*([A(i,1) A(i,2) A(i,3)]-P_0c)')';
    p.panels(f).P3=P_0c+(T*([A(i+M+1,1) A(i+M+1,2) A(i+M+1,3)]-P_0c)')';
    p.panels(f).P2=P_0c+(T*([A(i+M+2,1) A(i+M+2,2) A(i+M+2,3)]-P_0c)')';
    p.panels(f).P1=P_0c+(T*([A(i+1,1) A(i+1,2) A(i+1,3)]-P_0c)')';
    % Cslcolo le velocità nei vari punti
    Xb=(p.panels(f).P3(1)-p.panels(f).P4(1))/4+p.panels(f).P4(1);
    p.panels(f).A=P_0c+(T*([1000+Xb p.panels(f).P4(2) p.panels(f).P4(3)]-P_0c)')';
    p.panels(f).B=P_0c+(T*([Xb p.panels(f).P4(2) p.panels(f).P4(3)]-P_0c)')';
    Xc=(p.panels(f).P2(1)-p.panels(f).P1(1))/4+p.panels(f).P1(1);
    p.panels(f).D=P_0c+(T*([1000+Xc p.panels(f).P1(2) p.panels(f).P1(3)]-P_0c)')';
    p.panels(f).C=P_0c+(T*([Xc p.panels(f).P1(2) p.panels(f).P1(3)]-P_0c)')';
    %Normale
    p.panels(f).n=P_0c+(T*([0 -sin(d) cos(d)]-P_0c)')';
    %Punto di Controllo P
    l=sqrt((p.panels(f).P1(3)-p.panels(f).P4(3))^2+...
        (p.panels(f).P1(2)-p.panels(f).P4(2))^2)/2;
    Xb1=(p.panels(f).P3(1)-p.panels(f).P4(1))*3/4+p.panels(f).P4(1);
    Xc1=(p.panels(f).P2(1)-p.panels(f).P1(1))*3/4+p.panels(f).P1(1);
    tan_gamma=(Xc1-Xb1)/2;
    Z_pc=p.panels(f).P4(3)+l*sin(d);
    Y_pc=p.panels(f).P4(2)+b/4/M;
    X_pc=tan_gamma+Xb1;
    p.panels(f).P=P_0c+(T*([X_pc Y_pc Z_pc]-P_0c)')';
    f=f+1;
    i=i+1;
end
f=f-1;
%ora devo calcolare l'altra semiala; per farlo basta aggiungere gli stessi
%pannelli con tutte le y invertite, ed invertire i vari punti, poichè devo
%sempre usare le STESSE CONVENZIONI.
for i=1:f
    p.panels(f+i).P4(1)=p.panels(i).P1(1);
    p.panels(f+i).P4(2)=-p.panels(i).P1(2);
    p.panels(f+i).P4(3)=p.panels(i).P1(3);
    p.panels(f+i).P3(1)=p.panels(i).P2(1);
    p.panels(f+i).P3(2)=-p.panels(i).P2(2);
    p.panels(f+i).P3(3)=p.panels(i).P2(3);
    p.panels(f+i).P2(1)=p.panels(i).P3(1);
    p.panels(f+i).P2(2)=-p.panels(i).P3(2);
    p.panels(f+i).P2(3)=p.panels(i).P3(3);
    p.panels(f+i).P1(1)=p.panels(i).P4(1);
    p.panels(f+i).P1(2)=-p.panels(i).P4(2);
    p.panels(f+i).P1(3)=p.panels(i).P4(3);
    p.panels(f+i).A(1)=p.panels(i).D(1);
    p.panels(f+i).A(2)=-p.panels(i).D(2);
    p.panels(f+i).A(3)=p.panels(i).D(3);
    p.panels(f+i).B(1)=p.panels(i).C(1);
    p.panels(f+i).B(2)=-p.panels(i).C(2);
    p.panels(f+i).B(3)=p.panels(i).C(3);
    p.panels(f+i).C(1)=p.panels(i).B(1);
    p.panels(f+i).C(2)=-p.panels(i).B(2);
    p.panels(f+i).C(3)=p.panels(i).B(3);
    p.panels(f+i).D(1)=p.panels(i).A(1);
    p.panels(f+i).D(2)=-p.panels(i).A(2);
    p.panels(f+i).D(3)=p.panels(i).A(3);
    p.panels(f+i).P(1)=p.panels(i).P(1);
    p.panels(f+i).P(2)=-p.panels(i).P(2);
    p.panels(f+i).P(3)=p.panels(i).P(3);
    p.panels(f+i).n(1)=p.panels(i).n(1);
    p.panels(f+i).n(2)=-p.panels(i).n(2);
    p.panels(f+i).n(3)=p.panels(i).n(3);
end
f=f*2;
end
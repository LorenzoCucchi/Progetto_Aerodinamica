function [AIC]=AICMatrix(P)
%
% Input:
% P = struttura con i dati del profilo ottenuta con Panels
%
% Output:
% AIC = matrice A

N=length(P.d); % Numero di pannelli
AS=zeros(N);
bv=zeros(N,1);
bs=zeros(1,N);

bv_i_sum=0;
av=0;

for i=1:N
    n_i=P.n(i,:); % Normale del pannello i-esimo
    r_c_i=[P.x_c(i); P.y_c(i)]; % Posizione del centro del pannello i-esimo
    for j=1:N
        Q_j=[cos(P.theta(j)) sin(P.theta(j));
            -sin(P.theta(j)) cos(P.theta(j))];
        r_c_j=[P.x_c(j); P.y_c(j)];
        r_j_local=Q_j*(r_c_i-r_c_j);
        x_local=r_j_local(1);
        y_local=r_j_local(2);
        d_local=P.d(j);
        if j==i
            uj_ci_x_S_local=0; % Velocità autoindotta sul pannello j dalla sorgente del pannello j (componente x) sistema locale
            uj_ci_y_S_local=0.5; % Velocità autoindotta sul pannello j dalla sorgente del pannello j (componente y) sistema locale
        else
            uj_ci_x_S_local=1/(2*pi)*log((sqrt((x_local+d_local/2)^2+y_local^2))/(sqrt((x_local-d_local/2)^2+y_local^2))); % Velocità indotta sul pannello i dalla sorgente del pannello j (componente x) sistema locale
            uj_ci_y_S_local=1/(2*pi)*(atan((x_local+d_local/2)/y_local)-atan((x_local-d_local/2)/y_local)); % Velocità indotta sul pannello i dalla sorgente del pannello j (componente y) sistema locale
        end
        uj_ci_x_V_local=uj_ci_y_S_local; % Velocità indotta sul pannello i dal vortice del pannello j (componente x) sistema locale
        uj_ci_y_V_local=-uj_ci_x_S_local; % Velocità indotta sul pannello i dal vortice del pannello j (componente y) sistema locale
        uj_ci_S=Q_j'*[uj_ci_x_S_local; uj_ci_y_S_local]; % Velocità indotta sul pannello i dalla sorgente del pannello j 
        uj_ci_V=Q_j'*[uj_ci_x_V_local; uj_ci_y_V_local]; % Velocità indotta sul pannello i dal vortice del pannello j 
        AS(i,j)=dot(uj_ci_S,n_i); % Costruzione della Matrice A per le sorgenti
        bv_i_sum=bv_i_sum+dot(uj_ci_V,n_i); % Somma sulle j dei contributi di tutti i pannelli al termine del vortice
        if i==1 || i==N
            tau_i=P.tau(i,:); % Versore tangente al pannello i-esimo, calcolato solo per i=1 e i=N
            bs(j)=bs(j)+dot(uj_ci_S,tau_i); % Per i=1 e i=N calcola il nuovo contributo a bs(j) e lo somma
            av=av+dot(uj_ci_V,tau_i);
        end
    end
    bv(i)=bv_i_sum;
    bv_i_sum=0;
end

AIC=[AS bv;
    bs av];

end
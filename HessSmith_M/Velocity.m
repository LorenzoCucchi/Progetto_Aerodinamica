function [v]=Velocity(P,alpha,U,z)
%
% Input:
% P = struttura con i dati del profilo ottenuta con Panels
% alpha = AOA del profilo (deg)
% U = velocità asintotica
% z = soluzione del sistema lineare
%
% Output:
% v = velocità sul profilo

N=length(P.d);

U_inf=[U*cos(alpha) U*sin(alpha)];

U_S_sum=0;
U_V_sum=0;
v=[];

for i=1:N
    r_c_i=[P.x_c(i); P.y_c(i)]; % Posizione del centro del pannello i-esimo
    for j=1:N
        Q_j=[cos(P.theta(j)) sin(P.theta(j));
            -sin(P.theta(j)) cos(P.theta(j))]; % Matrice di rotazione del pannello j-esimo
        r_c_j=[P.x_c(j); P.y_c(j)]; % Posizione del centro del pannello j-esimo
        r_j_local=Q_j*(r_c_i-r_c_j); % r in coordinate locali
        x_local=r_j_local(1); % componente x di r in coordinate locali
        y_local=r_j_local(2); % componente y di r in coordinate locali
        d_local=P.d(j); % lunghezza del pannello j-esimo
        if j==i
            uj_ci_x_S_local=0; % Velocità autoindotta sul pannello j dalla sorgente del pannello j (componente x) sistema locale
            uj_ci_y_S_local=0.5; % Velocità autoindotta sul pannello j dalla sorgente del pannello j (componente y) sistema locale
        else
            uj_ci_x_S_local=1/(2*pi)*log((sqrt((x_local+d_local/2)^2+y_local^2))/(sqrt((x_local-d_local/2)^2+y_local^2))); % Velocità indotta sul pannello i dalla sorgente del pannello j (componente x) sistema locale
            uj_ci_y_S_local=1/(2*pi)*(atan((x_local+d_local/2)/y_local)-atan((x_local-d_local/2)/y_local)); % Velocità indotta sul pannello i dalla sorgente del pannello j (componente y) sistema locale
        end
        uj_ci_x_V_local=uj_ci_y_S_local; % Velocità indotta sul pannello i dal vortice del pannello j (componente x) sistema locale
        uj_ci_y_V_local=-uj_ci_x_S_local; % Velocità indotta sul pannello i dal vortice del pannello j (componente y) sistema locale
        uj_ci_S=(Q_j'*[uj_ci_x_S_local; uj_ci_y_S_local])'; % Velocità indotta sul pannello i dalla sorgente del pannello j 
        uj_ci_V=(Q_j'*[uj_ci_x_V_local; uj_ci_y_V_local])'; % Velocità indotta sul pannello i dal vortice del pannello j 
        U_S_sum=U_S_sum+z(j)*uj_ci_S;
        U_V_sum=U_V_sum+z(end)*uj_ci_V;
    end
    u_i=U_inf+U_S_sum+U_V_sum;
    v=[v; u_i];
    U_S_sum=0;
    U_V_sum=0;
end
    
function [d,gamma_scia]=Resistenza(rho,U_inf,Gamma_matrix_wing,Y,Z,M,p,alfa,beta,S)
%% calcolo resistenza
% HO una matrice Gamma_matrix di dimensione N*2M, calcolo i vortici in scia
% per cui prima sommo tutti i vortici rispetto ad N e ottengo il vortice a
% ferro di cavallo "Totale" visto dal bordo di uscita
Gamma_bordo_uscita=sum(Gamma_matrix_wing,1);
l_2d = rho*U_inf*Gamma_bordo_uscita;
Y_tot=[-fliplr(Y(1,2:end)) Y(1,:)];
Z_tot=[fliplr(Z(1,2:end)) Z(1,:)];
dy = abs(Y_tot(1,1)-Y_tot(1,2));

gamma_scia = [Gamma_bordo_uscita(1),-Gamma_bordo_uscita(1:(end-1))+Gamma_bordo_uscita(2:end),-Gamma_bordo_uscita(end)];
y = Y_tot(1,:);
z = Z_tot(1,:);

y_c = y(1:end-1)+0.5*dy;
z_c = 0.5*(z(1:end-1)+z(2:end));

v_ind = zeros(size(y_c));

for i = 1:length(y_c)
   for t = 1:length(gamma_scia)
   dist = sqrt((y_c(i)-y(t))^2+(z_c(i)-z(t))^2);        
   v_ind(i) = v_ind(i) + gamma_scia(t)*(-y_c(i)+y(t))/(4*pi*dist^2);
   end
end
y_vect=zeros(1,M);
for i=1:M
    y_vect(i)=p.panels(i).P(2);
end
y_vect=[-fliplr(y_vect) y_vect];
alfa_ind = -atan(v_ind/U_inf); % vettore incidenze indotte
figure ()
quiver(y_vect,zeros(1,2*M),zeros(1,2*M),v_ind,'color','b','linewidth',0.8), grid on
title (['Velocita indotta dalla scia ', '\alpha',' = ',...
    num2str(alfa),'°', '\beta',' = ', num2str(beta),'°'])
xlabel 'y', ylabel 'z'
d2d=l_2d.*sin(alfa_ind);
cd=sum(d2d)*dy/(0.5*rho*U_inf^2*S);
d=sum(d2d)*dy;
gamma_scia = [Gamma_bordo_uscita(1),-Gamma_bordo_uscita(1:(end-1))+Gamma_bordo_uscita(2:end),-Gamma_bordo_uscita(end)];

end
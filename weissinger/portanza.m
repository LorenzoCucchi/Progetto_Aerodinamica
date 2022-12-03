function [Cl_2D,L]=portanza(M,Gamma_matrix,rho,U_inf,c_med,d,b,N)
%fprintf("Calcolo portanza approccio 2D.... \n");
dy=b/2/M;
L=0;
L2d=zeros(2*M,1);
for f=1:2*M
    for j=1:N
        %Breve discorso sulle dimensioni: Gamma(Circolazione) è m^2/s, ho
        %una L2dF pari infatti a N/m, per cui per trovare il cl_2d devo
        %dividere per la corda media.
        L2d(f)=L2d(f)+rho*U_inf*Gamma_matrix(j,f)*cos(d);
    end
    Cl_2D(f)=L2d(f)/(0.5*rho*U_inf^2*c_med);
    L=L+L2d(f)*dy;
end
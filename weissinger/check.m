function check(p,FX,FY,FZ,Gamma,U)
% Funzione controllo per verficare che il risultato sia coerente con le
% condizioni al contorno.
fprintf("Controllo sulla soluzione.... \n");
[n,m]=size(Gamma);
V2=zeros(n,3);
V=[FX*Gamma FY*Gamma FZ*Gamma]+[U(1).*ones(n,1) U(2).*ones(n,1) U(3).*ones(n,1)];
for i=1:n
    % Le velocità sono corrette
    [V_tot]=Velocity(p,p.panels(i).P,Gamma)+U';
    V2(i,:)=V_tot;
    if dot(V_tot,p.panels(i).n)>=10^-12
        error("Velocità normale al profilo non trascurabile. PANNELLO n°%d",i);
    end
end
        
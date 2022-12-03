function [F,Moment]=force_2sup(p,Gamma,rho,U,G,M_tail,M)
[f,~]=size(Gamma);
% Ora calcolo il tutto rispetto alle due superfici 
Gamma_wing=Gamma(
F=[0 0 0];
Moment=[0 0 0];
for i=1:f_wing
    %%
    G_i=Gamma(i);
    %trovo la direzione di gamma, il contributo dei vortici a ferro di
    %cavallo ai lati del pannello è molto più piccolo rispetto agli altri,
    %ma per completezze ne tengo conto
    r_i=p.panels(i).C-p.panels(i).B;
    Middle_Point=0.5*(p.panels(i).C+p.panels(i).B);
    V=Velocity(p,Middle_Point,Gamma);
    F2=rho*G_i*(cross(V+U',r_i));
    F=F2+F;
    Moment=Moment+cross(Middle_Point-G,F2);
end
end
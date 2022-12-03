function [V_tot]=Velocity(p,P_vect,Gamma)
% Funzione che calcola la velocità indotta sul punto P_vect dai vortici a
% ferro di cavallo 
% INPUT:
% - p:    Struct contente tutte le info sui pannelli
% - P_vect:    Posizione punto 2
% - gamma:     Vettore contente le circolazione sui pannelli
% OUTPUT 
% - V:     velocità indotta
V_tot=[0 0 0];
[f,t]=size(Gamma);
for j=1:f
    I=Gamma(j);
    pA_j=p.panels(j).A;
    pB_j=p.panels(j).B;
    pC_j=p.panels(j).C;
    pD_j=p.panels(j).D;
    [V1]=BiotSavar(pA_j,pB_j,P_vect,I);
    [V2]=BiotSavar(pB_j,pC_j,P_vect,I);
    [V3]=BiotSavar(pC_j,pD_j,P_vect,I);
    V_tot=V1+V2+V3+V_tot;
end
end
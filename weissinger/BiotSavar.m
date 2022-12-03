function [V]=BiotSavar(p1,p2,c,g)
% Function per il calcolo della velocità indotta da un filamento vorticoso
% sul punto c. ATTENZIONE: gamma va da p1 a p2.
% INPUT:
% - p1:    Posizione punto 1
% - p2:    Posizione punto 2
% - c:     Posizione punto c
% OUTPUT 
% - V:     velocità indotta
% per evitare possibili "Nan" impongo che se la velocità si trova sul
% vortice allora la velocità indotta è nulla
r0=p2-p1;
r1=c-p1;
r2=c-p2;
distanza=norm(cross(r1,r2))/norm(r0);
toll=10^(-6);
if distanza < toll
  V= 0*r0;
else
  V=((g/(4*pi))*cross(r1,r2)/(norm(cross(r1,r2)))^2)*dot(r0,((r1/norm(r1))-(r2/norm(r2))));
end 
%R_cross=cross(r1,r2);
%V=g/4/pi*dot(r0,(r1/norm(r1)-r2/norm(r2)))*R_cross/norm(R_cross)^2;
end
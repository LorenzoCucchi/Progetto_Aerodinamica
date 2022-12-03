function GraficoIniziale(p,Gamma,U,d)
% Function "GraficoIniziale": Mostra la disposizione dei pannelli, la
% numerazione, le normali e le velocità tangenti. Inoltre verifica che le
% velocità e le normali siano perpendicolari.
% INPUT:
% - p:     Struct con le info dei pannelli
% - Gamma: Vettore contenete le circolazioni di ciascun vortice a ferro di
%          cavallo
% - U:     Velocità del flusso
% - d:     Se d==1 vengono numerati i pannelli
fprintf("Generazione Figura\n");
[f,t]=size(Gamma);
for i=1:f
    [V_tot]=Velocity(p,p.panels(i).P,Gamma)+U';
    if dot(V_tot,p.panels(i).n)>=10^(-10)
        error("ERROR: velocità calcolata sul pannelo %d viola le B.C.",i);
    end
    str=cellstr(num2str(i));
    figure (1)
    XP=[p.panels(i).P1(1) p.panels(i).P2(1) p.panels(i).P3(1) p.panels(i).P4(1) p.panels(i).P1(1)];
    YP=[p.panels(i).P1(2) p.panels(i).P2(2) p.panels(i).P3(2) p.panels(i).P4(2) p.panels(i).P1(2)];
    ZP=[p.panels(i).P1(3) p.panels(i).P2(3) p.panels(i).P3(3) p.panels(i).P4(3) p.panels(i).P1(3)];
    %Punto di Controllo
    plot3(p.panels(i).P(1),p.panels(i).P(2),p.panels(i).P(3),'*','color','b');
    hold on
    %plot3(p.panels(i).A(1)/100,p.panels(i).A(2),p.panels(i).A(3),'+','color','r');
    plot3(p.panels(i).B(1),p.panels(i).B(2),p.panels(i).B(3),'+','color','r');
    plot3(p.panels(i).C(1),p.panels(i).C(2),p.panels(i).C(3),'+','color','r');
    %Plot della normale
    quiver3(p.panels(i).P(1),p.panels(i).P(2),p.panels(i).P(3),p.panels(i).n(1),p.panels(i).n(2),p.panels(i).n(3),0.5,'color','g');
    quiver3(p.panels(i).P(1),p.panels(i).P(2),p.panels(i).P(3),V_tot(1),V_tot(2),V_tot(3),'color','b');
    %MarkerType
    %PLot del pannello
    plot3(XP,YP,ZP,'linewidth',0.1,'color','k');
    hold on
    if d==1
    text(p.panels(i).P(1),p.panels(i).P(2),p.panels(i).P(3),str)
    end
    
end
xlabel("X-axis");
ylabel("Y-axis");
zlabel("Z-axis");
grid on
axis equal
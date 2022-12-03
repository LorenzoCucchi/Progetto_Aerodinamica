function [Gamma,A,b,FX,FY,FZ]=LinearSystem(p,f,U)
% Aij=qij ed è la velocità indotta dal pannelo i dal j esimo ferro di
% cavallo, per cui ho un punto di controllo i esimo su cui calcolo la
% velocità indotta dal pannello j esimo, moltiplicata per la normale i
% esima del punto di controllo i esimo
A=zeros(f,f);
b=zeros(f,1);
for i=1:f
    %Definisco il centro c_i e la normale
    pP_i=p.panels(i).P;
    n_i=p.panels(i).n;
    for j=1:f
        pA_j=p.panels(j).A;
        pB_j=p.panels(j).B;
        pC_j=p.panels(j).C;
        pD_j=p.panels(j).D;
        [V1]=BiotSavar(pA_j,pB_j,pP_i,1);
        [V2]=BiotSavar(pB_j,pC_j,pP_i,1);
        [V3]=BiotSavar(pC_j,pD_j,pP_i,1);
        V_tot=V1+V2+V3;
        %F(i,j*3-2:j*3)=V_tot;
        FX(i,j)=V_tot(1);
        FY(i,j)=V_tot(2);
        FZ(i,j)=V_tot(3);
        A(i,j)=dot(V_tot,n_i);
    end
    b(i)=-dot(U,n_i);
end
Gamma=A\b;
check(p,FX,FY,FZ,Gamma,U);
end

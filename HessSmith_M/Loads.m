function [Cl,Cd,Cm]=Loads(P,Cp,U,alpha)
%
% Input:
% P = struttura con i dati del profilo ottenuta con Panels
% Cp = vettore con il Cp sui pannelli ottenuto da PressureCoeff
% U = velocità asintotica
% alpha = AOA del profilo (deg)
%
% Output:
% Cl, Cd, Cm = coefficienti di L, D e M del profilo

% Primo metodo per il calcolo di Cl: senza passare dal Cp
[A]=AICMatrix(P);
[b]=RHS(P,alpha,U);
x=A\b;
gamma=sum(P.d)*x(end);
Cl_1=2*gamma/U;

% Secondo metodo per il calcolo dei coefficienti: tramite il Cp
N=length(P.d);

U_vect = [U*cos(alpha); U*sin(alpha); 0];
alpha_vect = U_vect/U;
z = [0; 0; 1];
alpha_p = cross(z,alpha_vect);

Cl_2=0;
Cm=0;
Cd = 0;

for i=1:N
    n_i=[P.n(i,1) P.n(i,2) 0];
    r_c_i=[P.x_c(i)-1/4 P.y_c(i) 0];
    Cl_2=Cl_2-Cp(i)*P.d(i)*dot(n_i,alpha_p);
    Cd = Cd - Cp(i)*P.d(i)*dot(n_i,alpha_vect);
    Cm=Cm-Cp(i)*P.d(i)*dot(cross(r_c_i,n_i),z);
end

Cl=(Cl_1+Cl_2)/2;


end
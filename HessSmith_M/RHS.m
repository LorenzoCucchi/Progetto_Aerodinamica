function [b]=RHS(P,alpha,U)
%
% Input:
% P = struttura con i dati del profilo ottenuta con Panels
% alpha = AOA del profilo (deg)
% U = velocità asintotica
%
% Output:
% b = termine noto dell'equazione

N=length(P.d);

U_inf=[U*cos(alpha); U*sin(alpha)];

b=zeros(N+1,1);

for i=1:N
    b(i)=-dot(U_inf,P.n(i,:));
end

b(N+1)=-dot(U_inf,(P.tau(1,:)+P.tau(N,:)));

end
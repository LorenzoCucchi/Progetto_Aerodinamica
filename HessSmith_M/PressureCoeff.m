function [Cp]=PressureCoeff(v,U)
%
% Input:
% v = velocità sul profilo da Velocity
% U = velocità asintotica
%
% Output:
% Cp = Coefficiente di pressione sui pannelli

[N,~]=size(v);
Cp=[];

for i=1:N
    Cp=[Cp; 1-(norm(v(i,:))/U)^2];
end
end
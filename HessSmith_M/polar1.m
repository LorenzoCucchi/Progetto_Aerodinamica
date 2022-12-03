function [Cl,Cm,Cd,alpha_vec] = polar(P,U)


[A]=AICMatrix(P);
Cl = zeros(12,1);
Cd = zeros(12,1);
Cm = zeros(12,1);
alpha_vec = (-5:1:12)*pi/180;
for i = 1:length(alpha_vec)
    [b]=RHS(P,alpha_vec(i),U);
    z=A\b;
    [v]=Velocity(P,alpha_vec(i),U,z);
    [Cp]=PressureCoeff(v,U);
    [Cl(i),Cd(i),Cm(i)]=Loads(P,Cp,U,alpha_vec(i));
end


end

function A0=Calcolo_Alfa0(m,p)
%Definizione della linea media:
y_ml=@(s) m./p.^2.*(2.*p.*s-s.^2).*(s<p)+m./(1-p).^2.*((1-2.*p)+2.*p.*s-s.^2).*(s>=p);
yp_ml=@(s) 2.*m/((p)^2)*(p-s-1/2).*(s<p-1/2)+2.*m/((1-p).^2)*(p-s-1/2).*(s>=p-1/2);
% Definizione funzioni secondo la teoria dei profili sottili
Fpeso_1=@(s) sqrt((1/2+s)./(1/2-s));
dAlpha0=@(s) 2./pi.*Fpeso_1(s).*yp_ml(s);
% calcolo dei coefficienti
A0=integral(dAlpha0,-1/2,1/2);
end
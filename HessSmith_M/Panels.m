function [P]=Panels(x_p,y_p)
%
% Input:
% x_p, y_p = estremi dei pannelli così come escono da AirfoilShape
%
% Output:
% P = struttura con dati utili:
%     x_v = coordinate x degli estremi in ordine orario
%     y_v = coordinate y degli estremi in ordine orario
%     x_c = coordinate x dei centri in ordine orario
%     y_c = coordinate y dei centri in ordine orario
%     d = lunghezze dei pannelli
%     theta = orientamento dei pannelli
%     tau = versori tangenti
%     n = versori normali

P.x_v=[fliplr((x_p(:,2))') (x_p(2:end,1))']';
P.y_v=[fliplr((y_p(:,2))') (y_p(2:end,1))']';
P.x_c=(P.x_v(2:end)+P.x_v(1:end-1))/2;
P.y_c=(P.y_v(2:end)+P.y_v(1:end-1))/2;
P.d=sqrt((P.x_v(2:end)-P.x_v(1:end-1)).^2+(P.y_v(2:end)-P.y_v(1:end-1)).^2);
P.theta=wrapTo2Pi(atan2((P.y_v(2:end)-P.y_v(1:end-1))./P.d(1:end),(P.x_v(2:end)-P.x_v(1:end-1))./P.d(1:end)));
P.tau=[cos(P.theta) sin(P.theta)];
P.n=[-sin(P.theta),cos(P.theta)];

end
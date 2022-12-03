function [x_p,y_p]=AirfoilShape(NACA,N)
%
% Input:
% NACA = codice del profilo, a 4 o a 5 cifre, dato come stringa
% N = numero di pannelli desiderato
%
% Output:
% x_p, y_p = coordinate degli estremi dei pannelli del profilo, upper nelle
% prime colonne, lower nelle seconde

n=N/2+1; % Dato il numero di pannelli si calcolano i punti da mettere su x
x_vect=1-0.5*(1+cos((([1:n]-1)*pi)/(n-1)));
x_u=[];
x_l=[];
y_u=[];
y_l=[];

if length(NACA)==4
    m=str2num(NACA(1))/100;
    p=str2num(NACA(2))/10;
    SS=str2num(NACA(3:4));
    for x=x_vect
        y_t=5/100*SS*(0.2969*sqrt(x)-0.126*x-0.3516*x^2+0.2843*x^3-0.1036*x^4);
        if m==0 && p==0
            x_u_i=x;
            x_l_i=x;
            y_u_i=y_t;
            y_l_i=-y_t;
        else
            if x<=p
                y_c=m/(p^2)*(2*p*x-x^2);
                der=2*m/(p^2)*(p-x);
                theta=atan2(sin(der),(cos(der)));
            else
                y_c=m/((1-p)^2)*(1-2*p+2*p*x-x^2);
                der=2*m/((1-p)^2)*(p-x);
                theta=atan2(sin(der),(cos(der)));
            end
            x_u_i=x-y_t*sin(theta);
            x_l_i=x+y_t*sin(theta);
            y_u_i=y_c+y_t*cos(theta);
            y_l_i=y_c-y_t*cos(theta);
        end
        x_u=[x_u; x_u_i];
        x_l=[x_l; x_l_i];
        y_u=[y_u; y_u_i];
        y_l=[y_l; y_l_i];
    end 
elseif length(NACA)==5
    ddd=str2num(NACA(1:3));
    SS=str2num(NACA(4:5));
    switch ddd
        case 210
            k=361.4;
            q=0.0580;
        case 220
            k=51.64;
            q=0.1260;
        case 230
            k=15.957;
            q=0.2025;
        case 240
            k=6.643;
            q=0.2900;
        case 250
            k=3.230;
            q=0.3910;
        otherwise
            error('Profilo a 5 cifre non supportato')
    end
    for x=x_vect
        y_t=5/100*SS*(0.2969*sqrt(x)-0.126*x-0.3516*x^2+0.2843*x^3-0.1036*x^4);
        if x<=q
            y_c=k/6*(x^3-3*q*x^2+q^2*(3-q)*x);
            der=k/6*(3*x^2-6*q*x+q^2*(3-q));
            theta=atan2(sin(der),(cos(der)));
        else
            y_c=k/6*q^3*(1-x);
            der=-k/6*q^3;
            theta=atan2(sin(der),(cos(der)));
        end
        x_u_i=x-y_t*sin(theta);
        x_l_i=x+y_t*sin(theta);
        y_u_i=y_c+y_t*cos(theta);
        y_l_i=y_c-y_t*cos(theta);
        x_u=[x_u; x_u_i];
        x_l=[x_l; x_l_i];
        y_u=[y_u; y_u_i];
        y_l=[y_l; y_l_i];
    end
end

x_p=[x_u, x_l];
y_p=[y_u, y_l];

end

    
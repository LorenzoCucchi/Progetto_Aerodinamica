function [Cl, L_cl, L_circ, Gamma, Cdi, Di] = single_attitude(par_wing,par_field,N_plot,flag_plot)



%mapping
z=@(th) -par_wing.b/2 * cos(th);

N=5; %grado del polinomio di chebicev interpolante (introduce errore per troncamento)

B = Chev_coeffs(par_wing,z,N); %coefficienti del polinomio in theta

% theta_v=linspace(0,pi,N_plot); %vettore di theta dove valutare Gamma
z_v=linspace(-par_wing.b/2,par_wing.b/2,N_plot); %vettore di z dove valutare Gamma
%Valutazione del polinomio Gamma nei punti theta_v

Gamma=zeros(1,N_plot);

% for k=1:N_plot
%     summation=0;
%     for n=1:N
%         summation=summation + ( B(n) * sin(n*theta_v(k)) );
%     end
%     Gamma(k)=(2* par_wing.b * par_field.Uinf * summation) * par_wing.b/2 * sin(theta_v(k));
% 
% end

for k=1:N_plot
    summation=0;
    for n=1:N
        summation=summation + ( B(n) * sin(n* acos( -2/par_wing.b * z_v(k) )) );
    end
    Gamma(k)=(2* par_wing.b * par_field.Uinf * summation) ;

end



% calcolo della portanza
Cl=-pi* par_wing.AR * B(1);
L_cl=0.5*par_field.rho* par_field.Uinf^2 * par_wing.S * Cl;

Gamma_int=trapz(z_v,Gamma);
L_circ=-Gamma_int*par_field.rho*par_field.Uinf;

% calcolo della resistenza indotta
delta=0;
for n=2:N
    delta=delta + n * (B(n)/B(1))^2;
end


Cdi=(Cl^2 /(pi * par_wing.AR)) * (1+delta);
Di= 0.5* par_field.rho * par_field.Uinf^2 * par_wing.S * Cdi;

if (flag_plot)
    figure()
    plot(z_v,Gamma)
    %  axis equal
    xlim([-par_wing.b/2,par_wing.b/2])
    title('Gamma distribution over the wing')
   
    
    
    % [la portanza calcolata dalla teoria col Cl e la portanza generata tramite
    %    la circolazione dovrebbero risualtare uguali!]
    figure()
    l=-Gamma*par_field.rho*par_field.Uinf;
    plot(z_v,l)
    %  axis equal
    xlim([-par_wing.b/2,par_wing.b/2])
    title('Lift distribution over the wing')
   
end

end
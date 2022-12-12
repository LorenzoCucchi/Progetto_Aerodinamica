function [Cl,L_cl,L_circ,alfa_v] = alpha_plots(par_wing,par_field,alfa_inf,alfa_sup)

%plots parameters
N_plot=500;


alfa_v=linspace(alfa_inf,alfa_sup,50);

%mapping
z=@(th) -par_wing.b/2 * cos(th);



N=150; %grado del polinomio di chebicev interpolante (introduce errore per troncamento)
theta_v=linspace(0,pi,N_plot); %vettore di theta dove valutare il polinomio interpolante

%inizializzo i vettore per i plot
Cl=nan(1,length(alfa_v));
L_cl=nan(1,length(alfa_v));
L_circ=nan(1,length(alfa_v));

%genero i vettori in funzione di un cambio di alfa_g
for j=1:length(alfa_v)

    par_wing.alfa_g = @(z) deg2rad( alfa_v(j) );

    B = Chev_coeffs(par_wing,z,N); %coefficienti del polinomio in theta
    
    %definizione del polinomio inerpolante (la funzione da interpolare era la
    %circolazione)
    Gamma=zeros(1,N_plot);
    
    for k=1:N_plot
        summation=0;
        for n=1:N
            summation=summation + ( B(n) * sin(n*theta_v(k)) );
        end
        Gamma(k)=(2* par_wing.b * par_field.Uinf * summation) * par_wing.b/2 * sin(theta_v(k));
    
    end
    
    
    Cl(j)=-pi* par_wing.b^2 /par_wing.S * B(1);
    L_cl(j)=0.5*par_field.rho* par_field.Uinf^2 * par_wing.S * Cl(j);
    
    Gamma_int=trapz(linspace(0,pi,N_plot),Gamma);
    L_circ(j)=-Gamma_int*par_field.rho*par_field.Uinf;
    
    % [la portanza calcolata dalla teoria col Cl e la portanza generata tramite
    %    la circolazione dovrebbero risualtare uguali!]
    

end




end
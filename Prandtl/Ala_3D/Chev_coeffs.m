function B = Chev_coeffs(par_wing,z,N)
% funzione che calcola i coefficienti B del polinomio di Chebicev che fitta
% la distribuzione di circolazione data dalla geometria dell'ala

b=par_wing.b;
c=par_wing.c;
cl_a=par_wing.cl_a;
alfa_g=par_wing.alfa_g;
alfa_0=par_wing.alfa_0;

%scrittura sistema lineare Ain * Bn = Fi
thi=@(N,i) pi/N *(i-0.5); %punti di Chebicev


A=nan(N);
F=nan(N,1);
for i=1:N
    for n=1:N
        A(i,n)=( -4*b*sin(n*thi(N,i)) / c(z(thi(N,i))) ) - ( n*cl_a* sin(n*thi(N,i))/sin(thi(N,i)) );
    end
    F(i)=cl_a*(  alfa_g(z(thi(N,i)))  -  alfa_0(z(thi(N,i))) );
end

%risoluzione sistema lineare
B=A\F;

end
function visual_circ(X,Y,Z,Gamma,M,N)
% Visual_circ permette di fare un grafico, della distrubuzione di
% circolazione sulla superficie che si sta considerando
% INPUT 
% - X:  Matrice contente le x dei punti di estremità di ogni
%       pannello (M+1)*(N+1)
% - Z:  Matrice contente le x dei punti di estremità di ogni
%       pannello (M+1)*(N+1)
% - Y:  Matrice contente le x dei punti di estremità di ogni
%       pannello (M+1)*(N+1)
G1 =reshape(Gamma(1:(M*N)),M,N)';
G2 =reshape(Gamma((M*N)+1:2*M*N),M,N)';
surf(X,Y,Z,G1),colorbar, axis equal
hold on
surf(X,-Y,Z,G2),colorbar, axis equal
axis equal
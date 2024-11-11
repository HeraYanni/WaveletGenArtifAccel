function [G] = PSD_stat(zeta,N,omega,Sa,Ts,dom,n_iter)

p = 0.5; %

N_X = zeros(N,1);
eta_X = zeros(N,1);

% Spread factor
    delta_X = sqrt(1-(1/(1-zeta^2))*(1-(2/pi)*atan(zeta/sqrt(1-zeta^2)))^2);


for i=1:N
    % Mean zero crossing rate
    N_X(i) = (Ts/(2*pi))*omega(i)*(-log(p))^-1;
               
    % Peak factor
    eta_X(i) = sqrt(2*log(2*N_X(i)*(1-exp((-delta_X^1.2)*sqrt(pi*log(2*N_X(i)))))));
end


G = zeros(N,1);
Sum_G = 0 ;

% Start from j=2 to initialize because of omega(j-1) in the equation
for j=2:N
    % PSD calculation
      G(j) = ((2*zeta)/(omega(j)*pi-4*zeta*omega(j-1)))*((Sa(j)^2/eta_X(j)^2)-dom*Sum_G);
      Sum_G = Sum_G + G(j);
end


  %Target spectrum matching
for jj=1:n_iter
  
lamda0 = zeros(N,1);
lamda1 = zeros(N,1);
lamda2 = zeros(N,1);

N_X_lamda = zeros(N,1);
eta_X_lamda = zeros(N,1);
delta_X_lamda = zeros(N,1);
Sa_stoch = zeros(1,N);


for j=1:N
    omegaj=omega(j);
    
        n = 0;
    [lamda0(j)] = SpectralMoments(n,zeta,omegaj,omega,G);
    
        n = 1;
    [lamda1(j)] = SpectralMoments(n,zeta,omegaj,omega,G);
    
        n = 2;
    [lamda2(j)] = SpectralMoments(n,zeta,omegaj,omega,G);
    
    % Spread factor
    delta_X_lamda(j) = sqrt(1-lamda1(j)^2/(lamda0(j)*lamda2(j)));

    % Mean zero crossing rate
    N_X_lamda(j) = (Ts/(2*pi))*sqrt(lamda2(j)/lamda0(j))*(-log(p))^-1;
               
    % Peak factor
    eta_X_lamda(j) = sqrt(2*log(2*N_X_lamda(j)*(1-exp((-delta_X_lamda(j)^1.2)*sqrt(pi*log(2*N_X_lamda(j)))))));  
    
    Sa_stoch(j) = eta_X_lamda(j).*(omegaj.^2).*sqrt(lamda0(j));
    G(j) = G(j).*(Sa(j)./Sa_stoch(j)).^2;
end 

end

end
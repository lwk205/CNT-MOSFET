function [Grl,Grd,Gru,Gnl,Gnd,Gnu,Gpl,Gpd,Gpu,grL,ginL] = recursealg3d(Np,Al,Ad,Au,Sigin,Sigout)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        recursealg.m
% function [Grl,Grd,Gru,Gnl,Gnd,Gnu,Gpl,Gpd,Gpu] = recursealg3d(Np,Al,Ad,Au,Sigin,Sigout)
% recursive algorithm to solve for the diagonal elements only of 
% the Non-equilibrium Green's function
% HANDLES MATRICES BY 3 DIAGONALS
% Grl,Grd,Gru = retarded Green's function
% Gnl,Gnd,Gnu = electron Green's function
% Gpl,Gpd,Gpu = hole Green's function
% Np = size of the matrices
% Al,Ad,Au = matrix of coefficients
% Sigin = matrix of in-scattering self-energies (diagonal)
% Sigout = matrix of out-scattering self-energies (diagonal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_Gp = 'no';
Al_cr = conj(Au);
Ad_cr = conj(Ad);                                   % Hermitean conjugate of the coefficient matrix
Au_cr = conj(Al);
grL = zeros(1,Np);                                  % initialize left-connected function
ginL = zeros(1,Np);                                 % initialize left-connected in-scattering function
gipL = zeros(1,Np);                                 % initialize left-connected out-scattering function
Grl = zeros(1,Np-1);
Grd = zeros(1,Np);                                  % initialize the Green's function
Gru = zeros(1,Np-1);
Gnl = zeros(1,Np-1);
Gnd = zeros(1,Np);                                  % initialize the electron coherence function
Gnu = zeros(1,Np-1);
Gpl = zeros(1,Np-1);
Gpd = zeros(1,Np);                                  % initialize the hole coherence function
Gpu = zeros(1,Np-1);
grL(1)=1/Ad(1);                                     % step 1

for q=2:Np                                          % obtain the left-connected function
    grL(q)=1/(Ad(q)-Al(q-1)*grL(q-1)*Au(q-1));
end

gaL = conj(grL);                                        % advanced left-connected function
Grd(Np)=grL(Np);                                    % step 2

for q=(Np-1):-1:1
    Grl(q)=-Grd(q+1)*Al(q)*grL(q);                  % obtain the sub-diagonal of the Green's function
    Gru(q)=-grL(q)*Au(q)*Grd(q+1);                  % obtain the super-diagonal of the Green's function
%    Gru(q) = Grl(q);                               % exact for symmetric matrix A
    Grd(q)=grL(q)-grL(q)*Au(q)*Grl(q);              % obtain the diagonal of the Green's function
end

Gal = conj(Gru);
Gad = conj(Grd);                                    % advanced Green's function
Gau = conj(Gal);
ginL(1)=grL(1)*Sigin(1)*gaL(1);                     % step 3  

for q=2:Np
    sla2 = Al(q-1)*ginL(q-1)*Au_cr(q-1);
    prom = Sigin(q) + sla2;
    ginL(q)=grL(q)*prom*gaL(q);  % left-connected in-scattering function
end

Gnd(Np)=ginL(Np);                                   % step 4
Gnd = real(Gnd);

for q=(Np-1):-1:1
    Gnl(q) = - Grd(q+1)*Al(q)*ginL(q) - Gnd(q+1)*Al_cr(q)*gaL(q);                 
    % obtain the sub-diagonal of the electron Green's function
    Gnd(q) = ginL(q) + grL(q)*Au(q)*Gnd(q+1)*Al_cr(q)*gaL(q) - ( ginL(q)*Au_cr(q)*Gal(q) + Gru(q)*Al(q)*ginL(q) );
end

Gnu = conj(Gnl);                                    % super diagonal of the Green's function

switch flag_Gp
    case 'yes'
gipL(1)=grL(1)*Sigout(1)*gaL(1);                    % step 3  

for q=2:Np
    sla2 = Al(q-1)*gipL(q-1)*Au_cr(q-1);
    prom = Sigout(q) + sla2;
    gipL(q)=grL(q)*prom*gaL(q);  % left-connected in-scattering function
end

Gpd(Np)=gipL(Np);                                   % step 4
Gpd = real(Gpd);

for q=(Np-1):-1:1
    Gpl(q) = - Grd(q+1)*Al(q)*gipL(q) - Gnd(q+1)*Al_cr(q)*gaL(q);                 
    % obtain the sub-diagonal of the electron Green's function
    Gpd(q) = gipL(q) + grL(q)*Au(q)*Gpd(q+1)*Al_cr(q)*gaL(q) - ( gipL(q)*Au_cr(q)*Gal(q) + Gru(q)*Al(q)*gipL(q) );
end

Gpu = conj(Gpl);                                    % super diagonal of the hole function
    case 'no'
Gpl = 1i*(Grl-Gal) - Gnl;
Gpd = real(1i*(Grd-Gad) - Gnd);                            % hole Green's function
Gpu = 1i*(Gru-Gau) - Gnu;
end

jnzer = find(Gnd<0);
Gnd(jnzer) = 0;
jpzer = find(Gpd<0);
Gpd(jpzer) = 0;

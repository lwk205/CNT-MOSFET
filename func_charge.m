function [NIG] = func_charge(ee,Emcnt,AUD,Vd_bias,Vpp,kBT,Egh,WKB_flag,scat_flag,cont_flag,chorcu,maxfcnt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        func_charge.m
% function [NIG] = func_charge(ee,Emcnt,AUD,Vd_bias,Vpp,kBT,Egh,WKB_flag,chorcu,alpha)
% calculates charge density from the potential distribution by the
% Non-equilibrium Green's function method
% Nspec = spectrum of charge density (vs. energy)
% ee = energy
% Emcnt = energies of the middle of the bandgap (Np)
% AUD = off-diagonal element of the Hamiltonian (Np-1)
% Vd_bias = bias, source drain voltage
% Vpp = orbital coupling energy
% kBT = Boltzmann constant times the temperture
% Egh = half energy of the band gap
% WKB_flag = choice whether to calculate current by WKB or rigirously
% scat_flag = type of NEGF solution
% cont_flag = type of Source/Drain contacts
% chorcu = flag for charge vs, current
% maxfcnt = maximum # of energy grid points in myquad/myquad2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% in order to pass the values for populations
global eecount ii_band band_rec En_rec Gn_rec Gp_rec
global En_ex Gn_ex Gp_ex
global myEE myTr myCount Trraspr

amplesize = maxfcnt + 100;

Np=length(Emcnt);                       % size of the coordinate array
Sigma = sparse(Np,1);                   % initialize self energy
Gamma = sparse(Np,1);                   % initialize broadening rate
Sigma_in = sparse(Np,1);                % initialize in-scattering energy
Sigma_out = sparse(Np,1);               % initialize out-scattering energy
switch chorcu
    case 'charge'
        eta = 1e-14*1i;                   % infinitesimal imaginary number
    case 'current'
        eta = 1e-14*1i;
end

ep=ee+eta;                              % complex energy, with a small imaginary part
f_1=1/(1+exp(ee/kBT));                  % Fermi function from source
f_2=1/(1+exp((ee-(-Vd_bias))/kBT));     % Fermi function from drain
alphas=ep-Emcnt(1);                     % diagonal term at source
alphad=ep-Emcnt(Np);                    % diagonal term at drain
razn = sign(ee-Emcnt);                  % criterion of electron vs. hole conduction

switch cont_flag
    case 'semi'                         % self energy semi-semi contact
        b2=-AUD(1);                     % coupling to the contact, source
        Sigma_s=(-(b2*b2'-alphas^2-Vpp^2)+sqrt((b2*b2'-alphas^2-Vpp^2)^2-4*alphas^2*Vpp^2))/(2*alphas);
        b2=-AUD(Np-1);                  % coupling to the contact, drain
        Sigma_d=(-(b2*b2'-alphad^2-Vpp^2)+sqrt((b2*b2'-alphad^2-Vpp^2)^2-4*alphad^2*Vpp^2))/(2*alphad);
end

Gamma_s=1i*(Sigma_s-Sigma_s');               % injection rate, source
Gamma_d=1i*(Sigma_d-Sigma_d');               % injection rate, drain
Sigma_in_s = Gamma_s.*f_1;                  % influx through source
Sigma_in_d = Gamma_d.*f_2;                  % influx through drain
Sigma_out_s = Gamma_s.*(1-f_1);             % outflux through source
Sigma_out_d = Gamma_d.*(1-f_2);             % outflux through drain
Sigma_in(1) = Sigma_in(1) + Sigma_in_s;                     % start the lead self energies
Sigma_in(Np) = Sigma_in(Np) + Sigma_in_d;                   %
Sigma_out(1) = Sigma_out(1) + Sigma_out_s;                  %
Sigma_out(Np) = Sigma_out(Np) + Sigma_out_d;                %
Gamma(1) = Gamma(1) + Gamma_s;                              %
Gamma(Np) = Gamma(Np) + Gamma_d;                            %
Sigma(1) = Sigma(1) + Sigma_s;                              %
Sigma(Np) = Sigma(Np) + Sigma_d;                            %
Sigma_in0 = Sigma_in;                                       %
Sigma_out0 = Sigma_out;                                     %
Gamma0 = Gamma;                                             %
Sigma0 = Sigma;                                             % finish the lead self energies

switch scat_flag
   
    case 'ballist'                          % ballistic transport
        spB_s=sparse(Np,1);                         % initialized RHS
        spB_s(1)=1;                                 % such as to pick the source block of Green's function
        spB_d=sparse(Np,1);                         % initialized RHS
        spB_d(Np)=1;                                % such as to pick the drain block of Green's function
        Hdiag=Emcnt;                                % diagonal of hamiltonian
        Hdiag(1)=Hdiag(1)+Sigma_s;                  % add self energy source
        Hdiag(Np)=Hdiag(Np)+Sigma_d;                % add self energy drain
        Adiag = ee - Hdiag;                         % diagonal for the equation for Green's function (No broadening)
        AUD_up = [0; AUD];                          % upper diagonal for the equation for Green's function
        AUD_lo = [AUD; 0];                          % lower diagonal for the equation for Green's function
        A3diag = [AUD_lo Adiag AUD_up];             % three diagonals for the equation for Green's function
        G_inv = speye(Np,Np);                       % initialize space for the equation for Green's function
        G_inv = spdiags(A3diag,-1:1,G_inv);         % form the equation for Green's function
        G_s=G_inv\spB_s;                            % the source block of Green's function
        G_d=G_inv\spB_d;                            % the drain block of Green's function
        % Gamma = -2Im(Sigma) used explicitly here
        LDOS_S=-abs(G_s).^2*imag(Sigma_s)*2;        % local density of states, source
        LDOS_D=-abs(G_d).^2*imag(Sigma_d)*2;        % local density of states, drain
        % states filled by source and drain, ballistic regime, electrons if
        % ee>Emcnt, holes if ee<Emcnt
        LDOS = LDOS_S + LDOS_D;                                 % total densty of states
        Nelec = (f_2*LDOS_D+f_1*LDOS_S);                        % electron density
        Nhole = ((1-f_2)*LDOS_D+(1-f_1)*LDOS_S);                % hole density
        Nspec=0.5*(razn+1).*Nelec-0.5*(-razn+1).*Nhole;         % charge density
        Gnd = Nelec.';
        Gpd = Nhole.';
   
    case 'negf'                                                 % Non-equilibrium Green's function transport
        Sigma_in = Sigma_in0;       % total influx
        Sigma_out = Sigma_out0;   % total outflux
        Gamma = Gamma;                   % total relaxation rate
        Sigma = Sigma0;                   % total self energy
        AdiagS = ee - Emcnt - Sigma;                 % diagonal of the equation matrix
        % the universal routine to solve for the Green's function via the recursive (RGF) algorithm
        [Grl,Grd,Gru,Gnl,Gnd,Gnu,Gpl,Gpd,Gpu,grL,ginL] = recursealg3d(Np,AUD,AdiagS,AUD,Sigma_in,Sigma_out);
        Nspecn=real(Gnd).';                         % electron density
        Nspecp=real(Gpd).';                         % hole density
        Nspec = 0.5*(razn+1).*Nspecn-0.5*(-razn+1).*Nspecp;         % charge density
        Gad = conj(Grd);
        Aspec = i*(Grd - Gad);
        [res,oshi] = myisequal(Aspec.',Nspecn+Nspecp,1e-6);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate current %%%%%%%%%%%%%%%%
%iposit = floor(Np/2);                                  % where to calculate current
iposit = 1;
jcurr = 1:(Np-1);

switch chorcu
    case 'charge'
        NIG = Nspec;                             % charge density
    case 'current'
        NIG=zeros(Np+2,1);                          % placeholder for current and conductance
        switch WKB_flag
            case 'ballist'
                Tr=real(Gamma_d*(i*(G_d(Np)-G_d(Np)')-G_d(Np)*Gamma_d*G_d(Np)'));
            case 'rigor'
                Trinc = -i*AUD(iposit).*(Gnl(iposit)-Gnu(iposit));
                Trmass = -i*AUD(jcurr).'.*(Gnl(jcurr)-Gnu(jcurr));
                myEE(myCount) = ee;
                myTr(myCount) = Trinc;
                if isempty(Trraspr)
                    Trraspr = zeros(amplesize,Np-1);
                end
                Trraspr(myCount,:) = Trmass(1,:);
                myCount = myCount + 1;
        end
        
        switch WKB_flag
            case 'ballist'
                NIG(1)=Tr*(f_1-f_2);            % current spectrum
                NIG(2)=Tr/cosh(ee/(2*kBT))^2;   % conductance spectrum, normalized by 4e^2/h
                NIG(3:end) = Nspec;
            case 'rigor'
                NIG(1) = Trinc;
                NIG(2) = 0;         %Trinc/(f_1-f_2)/cosh(ee/(2*kBT))^2; not really useful for incoherent case
                NIG(3:end) = Nspec;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pass populations up %%%%%%%%%%
eecount = eecount + 1;
band_rec(eecount) = ii_band;
En_rec(eecount) = ee;

inzer = find(Gnd<0);
Gnd(inzer) = 0;

ipzer = find(Gpd<0);
Gpd(ipzer) = 0;

Gn_rec(eecount,:) = real(Gnd(1,:));
Gp_rec(eecount,:) = real(Gpd(1,:));


function [Fn_bias,Ne_bias,Id,Gd] = charge(Em,Vd_bias,Nband,delta,model,twob_flag,kBT,N1D,a,Vpp,n,Lab,...
    WKB_flag,scat_flag,cont_flag,chorcu,maxfcnt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        charge.m
% function [Fn_bias,Ne_bias,Id,Gd] = charge(Em,Vd_bias,Nband,delta,model,twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,chorcu,alpha)
% finds the charge density in 1D by solution of transport equation
% works for (n,0) zigzag CNT
% Fn_bias = dummy Fermi level, for the Poisson solver
% Ne_bias = charge density
% Id = source to drain current
% Gd = source to drain conductance
% Em = mid-gap energy
% Vd_bias = bias, source drain voltage
% Nband = number of sub-bands
% delta = parameter re. workfunction and temperature
% model = transport model
% '1' is the semiclassical ballistic transport; '2' is the quantum transport via NEGF.
% twob_flag = '0' - Only electron transport; '1' - both electron and hole transport
% kBT = thermal energy = k_B*T
% N1D = density constant for nanotubes
% a = grid size in Z direction
% Vpp = orbital coupling energy
% n = index of the nanotube
% Lab = reciprocal lattice period parameter
% WKB_flag = choice whether to calculate current by WKB or rigirously
% scat_flag = type of NEGF solution
% cont_flag = type of Source/Drain contacts
% chorcu = flag for charge vs, current
% maxfcnt = maximum # of energy grid points in myquad/myquad2
% Comments.
% Since energy is measured in eV here,
% Current = 2(for spin)*2(for valley degeneracy)*(q^2/h) [integral] dE
% transmssion * (fL-fR)
% Conductance = 2(for spin)*2(for valley degeneracy)*(q^2/h) [integral] dE
% transmssion * FT
% thermal broadening functon is
% FT = 1/(4kBT)(in this file) * 1/cosh^2(energy/2kBT)(n 'func_charge.m')


%% in order to pass the values for populations
global eecount ii_band En_rec band_rec Gn_rec Gp_rec
global En_old band_old Gn_old Gp_old
global En_ex Gn_ex Gp_ex


amplesize = maxfcnt + 100;
univconst                   % Define universal constants


Nw=length(Em);
Fn_bias=zeros(Nw,1);        % initialize ficticious fermi level
Ne_bias=zeros(Nw,1);        % initialize chage density
Ne_sub=zeros(Nw,1);         % initialize density in subbands
Id=0;                       % initialize drain current
Gd=0;                       % intialize drain conductance


switch model
    case 1                      % semiclassical ballistic approach
        for ii_band=1:Nband
            deltaband=delta*(6*ii_band-3-(-1)^ii_band)/4;
            [Ec_peak,i_peak]=max(Em+deltaband*kBT);
            for i_node=1:Nw
                eta=real(sqrt(((Ec_peak-Em(i_node))/kBT)^2-deltaband^2));
                zetas=(0-Em(i_node))/kBT;
                zetad=(-Vd_bias-Em(i_node))/kBT;
                if i_node<=i_peak
                    Ne_sub(i_node)=integral(zetas,zetad,eta,deltaband);
                elseif i_node>i_peak
                    Ne_sub(i_node)=integral(zetad,zetas,eta,deltaband);
                end
            end
            Ne_sub=Ne_sub*N1D;
            Ne_bias=Ne_bias+Ne_sub;
        end
        
        for ii_band=1:Nband
            if ii_band==1
                I1=round(n*2/3);
                b2=2*Vpp*cos(I1*pi/n);
                dI=sign(2/3-I1/n);
            else
                I=I1+dI;
                b2=2*Vpp*cos(I*pi/n);
                dI=sign(-dI)*ii_band;
            end
            Egh=abs(b2+Vpp);
            Ec_top=max(Em+Egh);
            Id_sub=log(1+exp((0-Ec_top)/kBT))-log(1+exp((-Vd_bias-Ec_top)/kBT));
            Id=Id+Id_sub;
            if Ec_top<0
                Gd=Gd+1;
            end
        end
        
        Id=(4*q^2*kBT/(2*pi*h_bar))*Id;
        
    case 2                     % NEGF
        Xe=2*Lab+a*[0:(Nw-1)];          % position of elements
        XI=[0 Xe Xe(Nw)+2*Lab];         % position of elements and boundary nodes
        Lw=a*(Nw-1)+4*Lab;
        EmXi=[Em(1); Em; Em(Nw)];
        % calculate the electrostatic potential matrix part of the Hamiltonian.
        Nab=floor(Lw/Lab);
        Np=2*Nab;
        Xcnt=zeros(Np,1);
        Em_atom=zeros(Np,1);
        
        for ii=1:Np
            Xcnt(ii)=acc*(3*ii/4+(1-(-1)^ii)/8-0.5);    % set atomistic grid inside coarse grid for Poisson
        end
        
        Em_atom=interp1(XI,EmXi,Xcnt);
        En=Em_atom;
        IG_tem=zeros(1,2);
        
        % initialize for population storage
        eecount = 0;                                    % call counter
        En_rec = zeros(amplesize,1);                    % placeholder for energies
        band_rec = zeros(amplesize,1);                  % placeholder for band indices
        Gn_rec = zeros(amplesize,Np);                   % placeholder for electron populations
        Gp_rec = zeros(amplesize,Np);                   % placeholder for electron populations
        % end initialization for population storage
        
        for ii_band=1:Nband
            if ii_band==1
                I1=round(n*2/3);
                b2=2*Vpp*cos(I1*pi/n);
                dI=sign(2/3-I1/n);
            else
                if mod(n,3)==0          % metallic
                    I=round(n*2/3)+ii_band-1;
                else                    % semiconducting
                    I=I1+dI;
                    dI=sign(-dI)*ii_band;
                end
                b2=2*Vpp*cos(I*pi/n);
            end
            
            Egh=abs(b2+Vpp);
            
            [E,E_least,E_most,E_peak,E_bot] = energygrid(twob_flag,Em_atom,Egh,Vd_bias);   % routine to set up energy grid
            % define the Hamiltonian
            Hup=zeros(Np-1,1);
            
            for ii=1:Np/2
                Hup(ii*2-1)=b2;
                if ii<Np/2
                    Hup(ii*2)=Vpp;
                end
            end
            
            AUD=-Hup;
            [En_ex,Gn_ex,Gp_ex] = sepBand(ii_band,En_old,band_old,Gn_old,Gp_old);

            chatol = 1e-8;  % integration accuracy for charge calculation , see matlab quad.m function
            curtol = 1e-8;  % integration accuracy for current calculation , see matlab quad.m function
            
            switch chorcu
                case 'charge'
                    [Ne_atom0 dum]=myquad(@func_charge,E_least,E_most,chatol,[],maxfcnt,Em_atom,AUD,Vd_bias,Vpp,kBT,Egh,...
                        WKB_flag,scat_flag,cont_flag,chorcu,maxfcnt);
                    Ne_atom=(4/(2*pi*(3*acc/4)))*Ne_atom0;
                    X_ring=0.5*(Xcnt(1:2:Np)+Xcnt(2:2:Np));
                    Ne_ring=0.5*(Ne_atom(1:2:Np)+Ne_atom(2:2:Np));
                    Ne_sub=interp1(X_ring,Ne_ring,Xe)';
                    Ne_bias=Ne_bias+Ne_sub;
                case 'current'
                    [IG_tem dum]=myquad2(@func_charge,E_least,E_most,curtol,chatol,[],maxfcnt,Em_atom,AUD,Vd_bias,Vpp,kBT,Egh,...
                        WKB_flag,scat_flag,cont_flag,chorcu,maxfcnt);
                    Id=Id+4*q^2/(2*pi*h_bar)*IG_tem(1);
                    Gd=Gd+(1/(4*kBT))*IG_tem(2);        % normalized by 4e^2/h
                    Ne_atom=(4/(2*pi*(3*acc/4)))*IG_tem(3:end);
                    X_ring=0.5*(Xcnt(1:2:Np)+Xcnt(2:2:Np));
                    Ne_ring=0.5*(Ne_atom(1:2:Np)+Ne_atom(2:2:Np));
                    Ne_sub=interp1(X_ring,Ne_ring,Xe)';
                    Ne_bias=Ne_bias+Ne_sub;
            end
            
            [En_rec,Gn_rec,Gp_rec] = sortEn(ii_band,En_rec,band_rec,Gn_rec,Gp_rec);         % sorting values of energy in a band
        end                                    %end of subbands
       
        switch chorcu
            case 'charge'
                eecount
                En_old = En_rec(1:eecount,:);
                band_old = band_rec(1:eecount,:);
                Gn_old = Gn_rec(1:eecount,:);   %% energy-position resolved electron distribution function
                Gp_old = Gp_rec(1:eecount,:);   %% energy-position resolved hole distribution function
            case 'current'
                count_curr = eecount
                En_old = En_rec(1:eecount,:);
                band_old = band_rec(1:eecount,:);
                Gn_old = Gn_rec(1:eecount,:);   %% energy-position resolved electron distribution function
                Gp_old = Gp_rec(1:eecount,:);   %% energy-position resolved hole distribution function
        end
end                                         % end of models
% convert charge density to dummy quasi-Fermi level for non-linear Poisson input

zeta=anti_dummy(Ne_bias'./N1D,delta,Nband); % 1st kind
Fn_bias=kBT*zeta+Em;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Title:      CNTDefineInputsShort.m
%   This file has the inputs for cylindrical CNT MOSFET simulation.
%   Parameters are defined at initialization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

univconst                   % Define universal constants
%% CNT related
n       = 10;
Lab     = 3/2*acc;  
d_cnt   = n*sqrt(3)*acc/pi
Vpp     = -3.0;     % tight-binding parameter
Egh1    = abs(2*Vpp*cos(round(n*2/3)*pi/n)+Vpp)
psi_md  = Egh1; % CNT work function
Te      = 400; %Lattice temperature
kBT     = k_B*Te/q;
delta   = abs(Egh1)/kBT;
N1D     = 8/(3*pi*acc*abs(Vpp))*kBT;

%Calculate Ne_0 and Ni, the charge density used in dummy and anti_dummy
Ne_0    = 4*kBT/(3*pi*abs(Vpp)*acc); %The charge density constants in nanotube carrier statistics
%end of calculating Ne_0 and Ni

%Gate contact material workfunction
phi_g   = Egh1; %Gate contact material workfuntion value
eps_ins = 16;   % oxide dielectric
t_ins   = 2e-9; % oxide thickness
Emgate  = phi_g - psi_md;

model   = 2;    %% '1' - for semiclassical ballistic transport; '2' - NEGF transport   
twob_flag= 1;   %% '0' - Only electron transport; '1' - both electron and hole transport
WKB_flag = 'ballist';   %% 'ballist' - NEGF formula for ballistic conduction  - use direct matrix inversion - fast for smaller devices
                        %% 'rigor' - NEGF for current calculation - use recursive Green's function (RGF) method - works for long devices (slow)
scat_flag = 'ballist';  %% Type of NEGF solution, 'ballist' - use direct matrix inversion; 
                        %% 'negf' - use recursive Green's function (RGF) method
                        
                        
%% NOTE: use either, (WKB_flag = 'rigor' && scat_flag = 'negf') combination, OR, (WKB_flag = 'ballist' && scat_flag = 'ballist') combination. 
cont_flag = 'semi';     %% S/D contact type for the self-energy, 'semi' - semi-conducting CNT contact;
maxfcnt = 10000; %% maximum number of energy grid point used in myquad.m and myquad2.m


%% NOTE: also see the 'chatol' and 'curtol' variables in charge.m for defining integration accuracy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nband   = 2;    %% Number of sub-bands for NEGF transport
a   = 2.5e-10;  %grid size in Z direction
b   = 1.0e-10;  %grid size in R direction
delta_x = 2*1e-9;    %% Source/Drain overlap with Gate
Lsd     = 3*1e-9;  %% Source/Drain length (without gate overlap)
Lg      = 10*1e-9;  %% Gate electrode length (includes source/drain overlap)
junction_l  = round((Lsd+delta_x)/a)+1;
junction_r  = round((Lsd+Lg-delta_x)/a)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define device boundaries for Poisson solution
col_nu      = round((2*Lsd+Lg)/a)+1;                    % number of points along the tube
row_nu      = round((d_cnt/2+t_ins)/b)+1;               % number of points in the radial direction
total_nu    = col_nu*row_nu;                            % total number of points
bound(1)    = 1;                                        % starting point
bound(2)    = col_nu;                                   % end point along tube
bound(3)    = round(d_cnt/2/b)*col_nu+1;                % 
bound(4)    = bound(3)+col_nu-1;
bound(5)    = round((d_cnt/2+t_ins)/b)*col_nu+1;
bound(6)    = bound(5)+round(Lsd/a);
bound(7)    = bound(6)+round(Lg/a);
bound(8)    = total_nu;
bound(9)    = junction_l;
bound(10)   = junction_r;
bound(11)   = col_nu;
bound(12)   = row_nu;
bound(13)   = total_nu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XI      = linspace(-(Lsd+Lg/2)*1e9,(Lsd+Lg/2)*1e9,col_nu); % in nanometers
% calculate model 2 related parameters
kf  = 2*pi/(3^1.5*acc);
t1  = 3*acc*abs(Vpp)/(4*a*sin(kf*a));

%% Doping densities
N_sd    = 1e9;   % source and drain diffusion 
N_body  = 0e9;      % body dopant concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Doping density array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nd  = zeros(bound(11),1);
for ii_node=1:round(bound(9)/2)
    Nd(ii_node)=N_sd*1;
end
for ii_node=round(bound(9)/2)+1:bound(9)
    Nd(ii_node)=N_sd;
end
for ii_node=bound(9)+1:bound(10)-1
    Nd(ii_node)=N_body;
end
for ii_node=bound(10):bound(10)+round(bound(9)/2)
    Nd(ii_node)=N_sd;
end
for ii_node=bound(10)+round(bound(9)/2)+1:bound(11)
    Nd(ii_node)=N_sd*1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Biasing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vg0     = -.5; %Gate contact bias
Vg_step = 0.1;
Ng_step = 12;  % number of steps
Vd0     = 0; %Drain contact bias
Vd_step = 0.05;
Nd_step = 12;  % number of steps


%% Self-consistency conditions
criterion_inner = 1e-7;                           %in eV
criterion_outer= 1e-3;                           %in eV
criterion_outer_semi = 3e-3;                      %in eV


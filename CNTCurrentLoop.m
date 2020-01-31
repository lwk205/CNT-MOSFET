function [Ie,Gd,Fn_bias,Ne_bias,Em_new,Em_cnt_new] = ...
    CNTCurrentLoop(Ne_bias_old,criterion_outer,Em_cnt_old,Vd,Nband,delta,model,...
    twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,scat_flag,cont_flag,...
    Em_old,criterion_inner,F_prime,Nd,eps_ins,bound,Egh1,b,Emgate,marquee,Err_File,maxfcnt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        CNTCurrentLoop.m
% function [Ie,Gd,Fn_bias,Ne_bias,Em,Em_cnt] = ...
%    CNTCurrentLoop(Ne_bias_old,criterion_outer,Em_cnt,Vd_bias(ii_vd),Nband,delta,model,...
%    twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,scat_flag,cont_flag,PhRel,PhRinel,homph,...
%    Em,criterion_inner,criterion_born,F_prime,Nd,eps_ins,bound,Egh1,b,Emgate,marquee,Err_File_H,alpha);
% Self-consistent current loop to determine current and charge after CNTChargeLoop.m has converged
% Ie = current
% Gd = conductance (for ballistic case)
% Fn_bias = fake fermi level for 'poisson'
% Ne_bias = new chagre densities
% Em_new = new mid-gap energy throught-out the device simulation grid points
% Em_cnt_new = new mid-gap energy at grid points on the CNT surface
% Em_old = old mid-gap energy throught-out the device simulation grid points
% Em_cnt_old = old mid-gap energy at grid points on the CNT surface
% Ne_bias_old = old charge distribution
% criterion_outer = tolerance for charge loop convergence
% Vd = source-drain bias
% Nband = number of bands
% delta = parameter re. workfunction and temperature
% model = NEGF or semiclassical
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
% criterion_inner = tolerance for the poisson loop convergence
% F_prime = jacobian for the poisson solver
% Nd = doping density array
% eps_ins =
% bound = array of integers for device boundaries
% Egh1 = energy gap
% b = grid size in radial direction
% Emgate = gate built-in potnetnial
% marquee = string - explanation for the error
% Err_File = handle of the error file
% maxfcnt = maximum # of energy grid points in myquad/myquad2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


iter_outer_lim = 14; %% Outer loop iteration limit
Ne_old = Ne_bias_old;
error_outer=1;
iter_outer=0;

fprintf(1,'Begin CNTCurrentLoop (NEGF)\n\n');
fprintf(Err_File,'Begin CNTCurrentLoop (NEGF)\n\n');

while((error_outer >= criterion_outer) & (iter_outer <= iter_outer_lim))
    [Fn_bias,Ne_bias,Ie,Gd] = charge(Em_cnt_old,Vd,Nband,delta,model,twob_flag,kBT,N1D,a,Vpp,n,Lab,...
        WKB_flag,scat_flag,cont_flag,'current',maxfcnt);

    [Em_new,Em_cnt_new,error_inner] = poisson(Fn_bias,Em_old,criterion_inner,F_prime,...
        Nd,eps_ins,kBT,Vpp,bound,Egh1,a,b,Emgate,N1D,delta,Nband);
    CNTerrWrite(Err_File,'Inner loop convergence criterion',error_inner);
    iter_outer=iter_outer+1;
    error_outer(iter_outer)=max(max(abs(Em_new-Em_old)))
    Em_old=Em_new;
    Em_cnt_old=Em_cnt_new;
end

if (iter_outer == iter_outer_lim + 1)
    fprintf(1,'CNTCurrentLoop reached limit: error_outer = %E , criterion_outer = %E\n',error_outer(iter_outer),criterion_outer);
    fprintf(Err_File,'CNTCurrentLoop reached limit: error_outer = %E , criterion_outer = %E\n',error_outer(iter_outer),criterion_outer);
end

CNTerrWrite(Err_File,marquee,error_outer);

fprintf(1,'CNTCurrentLoop (NEGF) finished\n\n');
fprintf(Err_File,'CNTCurrentLoop (NEGF) finished\n\n');



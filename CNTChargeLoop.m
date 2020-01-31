function [Fn_bias,Ne_bias,Em_new,Em_cnt_new] = ...
    CNTChargeLoop(criterion_outer,Em_cnt_old,Vd,Nband,delta,model,twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,scat_flag,cont_flag,...
                    Em_old,criterion_inner,F_prime,Nd,eps_ins,bound,Egh1,b,Emgate,marquee,Err_File,maxfcnt)
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        CNTChargeLoop.m 
% function [Fn_bias,Ne_bias,Em_new,Em_cnt_new,Em_old,Em_cnt_old] = ...
%    CNTChargeLoop(criterion_outer,Em_cnt_old,Vd,Nband,delta,model,twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,...
%    Em_old,criterion_inner,F_prime,Nd,eps_ins,bound,Egh1,b,Emgate,marquee,Err_File,alpha
% self-consistent loop to calculate the self-consistent solution for the 
% transport and electric field
% Fn_bias = fake fermi level for 'poisson'
% Ne_bias = chagre densities
% Em_new = new mid-gap energy throught-out the device simulation grid points 
% Em_cnt_new = new mid-gap energy at grid points on the CNT surface
% Em_old = old mid-gap energy throught-out the device simulation grid points 
% Em_cnt_old = old mid-gap energy at grid points on the CNT surface
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

iter_outer_lim = 5; %% Outer loop iteration limit
%%%%%%%%%%%%%%%% Ballistic resultsd %%%%%%%%%%%%%%%
error_outer=1;
iter_outer=0;

while(error_outer >= 10*criterion_outer & iter_outer <= iter_outer_lim)
    [Fn_bias,Ne_bias,Id,Gd] = charge(Em_cnt_old,Vd,Nband,delta,model,twob_flag,kBT,N1D,a,Vpp,n,Lab,...
        WKB_flag,scat_flag,cont_flag,'charge',maxfcnt);

    [Em_new,Em_cnt_new,error_inner] = poisson(Fn_bias,Em_old,criterion_inner,F_prime,...
        Nd,eps_ins,kBT,Vpp,bound,Egh1,a,b,Emgate,N1D,delta,Nband);
    CNTerrWrite(Err_File,'Inner loop convergence criterion',error_inner);
    iter_outer=iter_outer+1;
    error_outer(iter_outer)=max(max(abs(Em_new-Em_old)))
    Em_old=Em_new;
    Em_cnt_old=Em_cnt_new;
end

if ((iter_outer == iter_outer_lim + 1) & (model == 2))
    fprintf(1,'CNTChargeLoop reached limit: error_outer = %E , criterion_outer = %E\n',error_outer(iter_outer),criterion_outer);
    fprintf(Err_File,'CNTChargeLoop reached limit: error_outer = %E , criterion_outer = %E\n',error_outer(iter_outer),criterion_outer);
end

CNTerrWrite(Err_File,marquee,error_outer);

if (model == 2)
    fprintf(1,'Ballistic CNTChargeLoop (NEGF) finished\n\n');
    fprintf(Err_File,'Ballistic CNTChargeLoop (NEGF) finished\n\n');
end
%%%%%%%%%%%%%%% End Ballistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

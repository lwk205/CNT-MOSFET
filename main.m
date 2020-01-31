%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        main.m
% The top level program for carbon nanotube (MOSFET) transistor



clc;
clear all;
close all;

global En_old band_old Gn_old Gp_old;
global myEE myTr myCount Trraspr

CNTInitialWrite;                                    % write to file in the beginning
CNTDefineInputs;                                    % inputs, for regular simulation
F_prime = CNTDefineFprime(bound,a,b,eps_ins);       % Specify F_prime matrix for Jacobian in Poisson solver
CNTChargeNeut;                                      % charge neutrality for initial guess
for ii_vg = 1:Ng_step+1                             % start Vg loop
    Vg_bias(ii_vg) = Vg0 + Vg_step*(ii_vg-1);       % set gate voltage
    Emgate = -Vg_bias(ii_vg) + phi_g - psi_md;      % set energy midgap due to gate
    for ii_vd = 1:Nd_step+1                         % start Vds loop
        Vd_bias(ii_vd) = Vd0 + Vd_step*(ii_vd-1);   % set source-drain voltage
        mu_d = -Vd_bias(ii_vd);                     % set chemical potnetial at drain
        
        fprintf(Err_File_H,'\n**** Bias Point: Vg = %.3f V, Vd = %.3f V ****\n',Vg_bias(ii_vg),Vd_bias(ii_vd));
        fprintf(1,'\n**** Bias Point: Vg = %.3f V, Vd = %.3f V ****\n',Vg_bias(ii_vg),Vd_bias(ii_vd));
        Fn_bias = zeros(bound(11),1);               % initialize dummy fermi level for 'poisson'
        Em  = zeros(bound(13),1);               % initialize electric field
        [Em,Em_cnt,error_inner] = poisson(Fn_bias,Em,criterion_inner,F_prime,...
            Nd,eps_ins,kBT,Vpp,bound,Egh1,a,b,Emgate,N1D,delta,Nband);      % initial solution for electric field
        marquee = 'Semi-classical solution: Outer loop convergence criterion';
        
        [Fn_bias,Ne_bias,Em,Em_cnt]  = ...   % self-consistent loop: semiclassical as an initial guess
            CNTChargeLoop(criterion_outer_semi,Em_cnt,Vd_bias(ii_vd),Nband,delta,1,...
            twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,scat_flag,cont_flag,...
            Em,criterion_inner,F_prime,Nd,eps_ins,bound,Egh1,b,Emgate,marquee,Err_File_H,maxfcnt);
        marquee = 'CNTChargeLoop: Outer loop convergence criterion';
        
        [Fn_bias,Ne_bias,Em,Em_cnt]  = ...   % self-consistent loop: NEGF to determine charge
            CNTChargeLoop(criterion_outer,Em_cnt,Vd_bias(ii_vd),Nband,delta,model,...
            twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,scat_flag,cont_flag,...
            Em,criterion_inner,F_prime,Nd,eps_ins,bound,Egh1,b,Emgate,marquee,Err_File_H,maxfcnt);
        marquee = 'CNTCurrentLoop: Outer loop convergence criterion';
        
        [Ie_tem,Gd_tem,Fn_bias,Ne_bias,Em,Em_cnt] = ... % self-consistent current loop to determine current and charge
            CNTCurrentLoop(Ne_bias,criterion_outer,Em_cnt,Vd_bias(ii_vd),Nband,delta,model,...
            twob_flag,kBT,N1D,a,Vpp,n,Lab,WKB_flag,scat_flag,cont_flag,...
            Em,criterion_inner,F_prime,Nd,eps_ins,bound,Egh1,b,Emgate,marquee,Err_File_H,maxfcnt);

        Ne(:,ii_vg,ii_vd)=Ne_bias;
        Ec(:,ii_vg,ii_vd)=Em_cnt+Egh1;
        Ie(ii_vg,ii_vd)=Ie_tem;
        E_sub = Ec;
        Ne_sub = Ne;
        Trraspr = Trraspr(1:myCount-1,:);
        
        % save partial results (In case simulation stops)
%         save mos_results XI E_sub n Ne_sub Vd_bias Vg_bias phi_g Ie Ng_step Nd_step Lg t_ins Egh1 ...
%         En_old Gn_old Gp_old myEE myTr Trraspr I_continuity

        save mos_results XI E_sub n Ne_sub Vd_bias Vg_bias phi_g Ie Ng_step Nd_step Lg a t_ins Egh1 

        Ie
    end                                             % end Vds loop
end                                                 % end Vg loop

CNTSaveResults;                                     % save results to file

%CNTPlotResults; 


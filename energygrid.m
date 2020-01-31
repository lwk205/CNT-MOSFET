function [E,E_least,E_most,E_peak,E_bot] = energygrid(twob_flag,Em_atom,Egh,Vd_bias)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [E,E_least,E_most,E_peak] = energygrid(twob_flag,Em_atom,Egh,Vd_bias)
% Title:        energygrid.m 
% set up energy grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_step=1e-3;

if twob_flag==0
     Ef_tail_up=0.4;
     Ef_tail_low=0.2;
     E_peak=max(0,max(Em_atom+Egh));
     E_bot=min(Em_atom+Egh);
else
     Ef_tail_up=0.4;
     Ef_tail_low=0.2;
     E_peak=max(0,max(Em_atom+Egh));
     E_bot=min(-Vd_bias,min(Em_atom-Egh));
end

E_least = E_bot-Ef_tail_low;
E_most = E_peak+Ef_tail_up;

num_pred = round((E_most-E_least)/E_step);
E_number= num_pred + mod(num_pred,2) + 1;
E=linspace(E_least,E_most,E_number);
delta_E=(E_most-E_least)/(E_number-1);
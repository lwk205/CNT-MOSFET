function [F]=dummy(zeta,delta,der_flag,Nband)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        dummy.m
% function [F]=dummy(zeta,delta,der_flag,Nband)
% This function computes the carrier statistical integral for CNT square root E(k)
% input: zeta=(mu-Em)/kT; delta=Egh/kT; der_flag=0 for integral, =1 for dI/d(zeta)
% F =
% zeta = 
% delta =
% der_flag
% Nband

tail_up=15;
F=zeros(1,length(zeta));

for ii_band=1:Nband
    deltaband=delta*(6*ii_band-3-(-1)^ii_band)/4;
    lim_up=max(max(abs(zeta)),deltaband)+tail_up;
    
    if der_flag==0
%       ff=inline('1./(1+exp(sqrt(x^2+delta^2)-zeta))-1./(1+exp(sqrt(x^2+delta^2)+zeta))','x','zeta','delta');
%       F=F+myquad(ff,0,lim_up,1e-7,[],zeta,deltaband);
        F=F+myquad('dufu1',0,lim_up,1e-7,[],10000,zeta,deltaband);  %% use myquad with maxfcnt = 10000
    elseif der_flag==1
%        ff=inline('0.25./cosh((sqrt(x^2+delta^2)-zeta)./2).^2+0.25./cosh((sqrt(x^2+delta^2)+zeta)./2).^2','x','zeta','delta');
%        F=F+myquad(ff,0,lim_up,1e-7,[],zeta,deltaband);
        F=F+myquad('dufu2',0,lim_up,1e-7,[],10000,zeta,deltaband);   %% use myquad with maxfcnt = 10000
    end
end


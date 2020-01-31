% This function computes the carrier statistical integral for CNT square root E(k)
% input: zeta=(mu-Em)/kT; delta=Egh/kT; der_flag=0 for integral, =1 for dI/d(zeta)

function [F]=integral(zetas,zetad,eta,delta)
infs=1e-5;

tail_up=20;
dx0=1e-2;
lim_ups=max(zetas,delta)+tail_up;
dx=dx0; nx=round(lim_ups/dx)+2;
xx=linspace(0,lim_ups,nx); dx=xx(2)-xx(1);
yy=1./(1+exp(sqrt(xx.^2+delta^2)-zetas));
F1=(sum(yy)-yy(1)/2-yy(nx)/2)*dx;

clear xx yy nx
lim_upd=max(eta,zetad)+tail_up;
dx=dx0; nx=round((lim_upd-eta)/dx)+2;
xx=linspace(eta,lim_upd,nx); dx=xx(2)-xx(1);
yy=1./(1+exp(sqrt(xx.^2+delta^2)-zetad));
F3=(sum(yy)-yy(1)/2-yy(nx)/2)*dx;

if eta==0   % get rid of the warning message at the top of the barrier
    eta=infs;
end 

clear xx yy nx

dx=dx0; nx=round(eta/dx)+2;
xx=linspace(0,eta,nx); dx=xx(2)-xx(1);
yy=1./(1+exp(sqrt(xx.^2+delta^2)-zetas));
F2=(sum(yy)-yy(1)/2-yy(nx)/2)*dx;

F=(F1+F2+F3)*0.5;



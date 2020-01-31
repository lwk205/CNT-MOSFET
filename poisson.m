function [Em,Em_cnt,error_inner] = poisson( Fn,Em_old,criterion_inner,F_prime,...
                                            Nd,eps_ins,kBT,Vpp,bound,Egh1,a,b,Emgate,N1D,delta,Nband)
                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title:        poisson.m
% function [Em,Em_cnt,error_inner] = poisson(   Fn,Em_old,criterion_inner,F_prime,...
%                                               Nd,eps_ins,kBT,Vpp,bound,Egh1,a,b,Emgate,N1D,delta,Nband)
% This program solves  2-D Poisson equation in cylindrical coordinates using Newton-Ralphson method
% for the coaxially gated CNTFETs
% Em = new mid-gap potential through-out the device grid points
% Em_cnt = new mid-gap potential at grid points on the CNT surface.
% error_inner = 
% Fn = 
% Em_old = old mid-gap potential through-out the device grid points
% criterion_inner =
% F_prime = 
% Nd = 
% eps_ins = 
% kBT = 
% Vpp = 
% bound = 
% Eghl = 
% a = 
% b = 
% Emgate = 
% N1D = 
% delta = 
% Nband = 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


univconst                   % Define universal constants
delta_Em=zeros(bound(13),1);
F=zeros(bound(13),1);
MF_prime=sparse(bound(13),bound(13));
Em=zeros(bound(13),1);
dummy_fun=zeros(bound(11),1);
dummy_fun_prime=zeros(bound(11),1);
iter_inner=0;
error_inner=1;
Em=real(Em_old);
Em_cnt=Em(bound(3):bound(4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% STRAT OF INNER LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(error_inner>=criterion_inner)
iter_inner=iter_inner+1;
%dummy specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zetan=(Fn-Em_cnt)./kBT;
dummy_fun=N1D*dummy(zetan,delta,0,Nband);
dummy_fun_prime=-(N1D/kBT)*dummy(zetan,delta,1,Nband);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate F vector. (residue for Newton-Ralphson method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify the F vector in the air region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii_node=bound(2)+1: bound(3)-1
    %Calculate the discritized R position
    mm=floor((ii_node-1)/bound(11));
    F(ii_node)=-(a/b)*(1-1/(2*mm))*Em(ii_node-bound(11))...
        -(b/a)*Em(ii_node-1)...
        +2*(a/b+b/a)*Em(ii_node)...
        -(b/a)*Em(ii_node+1)...
        -(a/b)*(1+1/(2*mm))*Em(ii_node+bound(11));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify the F vector in the gate insulator region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii_node=bound(4)+1: bound(5)-1
    %Calculate the discritized R position
    mm=floor((ii_node-1)/bound(11));
    F(ii_node)=-(a/b)*eps_ins*(1-1/(2*mm))*Em(ii_node-bound(11))...
        -(b/a)*eps_ins*Em(ii_node-1)...
        +2*(a/b+b/a)*eps_ins*Em(ii_node)...
        -(b/a)*eps_ins*Em(ii_node+1)...
        -(a/b)*eps_ins*(1+1/(2*mm))*Em(ii_node+bound(11));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify the F vector on the nanotube shell (gate insulator/air interface).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii_node=bound(3):bound(4)
    %Calculate the discritized R position
    mm=floor((ii_node-1)/bound(11));
    nn=mod((ii_node-1),bound(11))+1;
    F(ii_node)=-(1/8)*(b/a)*(1-1/(4*mm))*Em(ii_node-bound(11)-1)...
        -((a/b)*(1-1/(2*mm))-(1/4)*(b/a)*(1-1/(4*mm)))*Em(ii_node-bound(11))...
        -(1/8)*(b/a)*(1-1/(4*mm))*Em(ii_node-bound(11)+1)...
        -(3/8)*(b/a)*((1-1/(4*mm))+eps_ins*(1+1/(4*mm)))*Em(ii_node-1)...
        +(eps_ins*(a/b)*(1+1/(2*mm))+(a/b)*(1-1/(2*mm))+(3/4)*(b/a)*eps_ins*(1+1/(4*mm))+(3/4)*(b/a)*(1-1/(4*mm)))*Em(ii_node)...
        -(3/8)*(b/a)*((1-1/(4*mm))+eps_ins*(1+1/(4*mm)))*Em(ii_node+1)...
        -(1/8)*(b/a)*(1+1/(4*mm))*eps_ins*Em(ii_node+bound(11)-1)...
        -((a/b)*(1+1/(2*mm))-(1/4)*(b/a)*(1+1/(4*mm)))*eps_ins*Em(ii_node+bound(11))...
        -(1/8)*(b/a)*(1+1/(4*mm))*eps_ins*Em(ii_node+bound(11)+1)...
        +(a*q/(2*pi*eps_o*mm*b))*(Nd(nn)-dummy_fun(nn));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify the F vector on boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The R=0 boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii_node=bound(1):bound(2)
    F(ii_node)=Em(ii_node)-Em(ii_node+bound(11));
end
%The gate contact
for ii_node=bound(5):bound(8)
    if ii_node<bound(6)
        F(ii_node)=Em(ii_node)-Em(ii_node-bound(11));
    elseif (ii_node>=bound(6))&(ii_node<=bound(7))
        F(ii_node)=Em(ii_node)-Emgate;
    elseif ii_node>bound(7)
        F(ii_node)=Em(ii_node)-Em(ii_node-bound(11));
    end
end

%The left and right boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii_node_l=bound(1);
ii_node_r=bound(2);

for jjj=1:bound(12)
    if jjj==1
        F(ii_node_l)=Em(ii_node_l)-(1/2)*Em(ii_node_l+1)-(1/2)*Em(ii_node_l+bound(11));
        F(ii_node_r)=Em(ii_node_r)-(1/2)*Em(ii_node_r-1)-(1/2)*Em(ii_node_r+bound(11));
    elseif (jjj>1)&(jjj<bound(12))
        F(ii_node_l)=Em(ii_node_l)-Em(ii_node_l+1);
        F(ii_node_r)=Em(ii_node_r)-Em(ii_node_r-1);
    elseif jjj==bound(12)
        F(ii_node_l)=Em(ii_node_l)-(1/2)*Em(ii_node_l+1)-(1/2)*Em(ii_node_l-bound(11));
        F(ii_node_r)=Em(ii_node_r)-(1/2)*Em(ii_node_r-1)-(1/2)*Em(ii_node_r-bound(11));
    end
    ii_node_l=ii_node_l+bound(11);
    ii_node_r=ii_node_r+bound(11);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MF_prime=sparse(bound(13),bound(13));

for ii_node=bound(3)+1:bound(4)-1
    nn=mod((ii_node-1),bound(11))+1;
    MF_prime(ii_node,ii_node)=-(a*q/(2*pi*eps_o*mm*b))*dummy_fun_prime(nn);
end

MF_prime=F_prime+sparse(MF_prime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF EVALUATING MF_prime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVING FOR delta_Ec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta_Ec=-sparse(MF_prime)\sparse(F);
delta_Ec=-MF_prime\F;       % matrices on the right are already sparse
dsel = delta_Ec(1:bound(13));
dselab = abs(dsel);
jf1 = find(1.0<dselab & dselab<3.7);
dsel(jf1)=sign(dsel(jf1)).*dselab(jf1).^(1/5);
jf2 = find(dselab>=3.7);
dsel(jf2)=sign(dsel(jf2)).*log(dselab(jf2));
delta_Ec(1:bound(13)) = dsel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SOLVING FOR delta_Ec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Em=Em+delta_Ec;
Em_cnt=Em(bound(3):bound(4));
error_inner(iter_inner)=max(abs(real(F)));
max_delta_Ec=max(abs(real(full(delta_Ec))));
MF_prime=sparse(bound(13),bound(13));
F=zeros(bound(13),1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End OF INNER LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
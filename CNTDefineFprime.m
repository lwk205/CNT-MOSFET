function    F_prime = CNTDefineFprime(bound,a,b,eps_ins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Title:      CNTDefineFprime.m
%   This defines the F_prime matrix used in finite-difference poisson
%   solver
%   function    F_prime = CNTDefineFprime(bound,a,b,eps_ins)
%   bound   = devise boundary values
%   a       = grid size in Z direction
%   b       = grid size in R direction
%   eps_ins = insulator dielectric const
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b13 = bound(13);
F_prime = sparse(b13,b13);
rowzero = sparse(1,b13);

%The air region
for ii_node=bound(2)+1:bound(3)-1
    mm=floor((ii_node-1)/bound(11));
    F_prime(ii_node,ii_node-bound(11))=-(a/b)*(1-1/(2*mm));
    F_prime(ii_node,ii_node-1)=-b/a;
    F_prime(ii_node,ii_node)=2*(a/b+b/a);
    F_prime(ii_node,ii_node+1)=-b/a;
    F_prime(ii_node,ii_node+bound(11))=-(a/b)*(1+1/(2*mm));
end

%The gate insulator region
for ii_node=bound(4)+1:bound(5)-1
    mm=floor((ii_node-1)/bound(11));
    F_prime(ii_node,ii_node-bound(11))=-(a/b)*(1-1/(2*mm))*eps_ins;
    F_prime(ii_node,ii_node-1)=-b/a*eps_ins;
    F_prime(ii_node,ii_node)=2*(a/b+b/a)*eps_ins;
    F_prime(ii_node,ii_node+1)=-b/a*eps_ins;
    F_prime(ii_node,ii_node+bound(11))=-(a/b)*(1+1/(2*mm))*eps_ins;
end

%The carbon nanotube shell (Air/Gate insulator interface)
for ii_node=bound(3):bound(4)
    mm=floor((ii_node-1)/bound(11));
    nn=mod((ii_node-1),bound(11))+1;
    F_prime(ii_node,ii_node-bound(11)-1)=-(1/8)*(b/a)*(1-1/(4*mm));
    F_prime(ii_node,ii_node-bound(11))=-((a/b)*(1-1/(2*mm))-(1/4)*(b/a)*(1-1/(4*mm)));
    F_prime(ii_node,ii_node+bound(11)+1)=-(1/8)*(b/a)*(1-1/(4*mm));
    F_prime(ii_node,ii_node-1)=-(3/8)*(b/a)*((1-1/(4*mm))+eps_ins*(1+1/(4*mm)));
    F_prime(ii_node,ii_node)=+(eps_ins*(a/b)*(1+1/(2*mm))+(a/b)*(1-1/(2*mm))+(3/4)*(b/a)*eps_ins*(1+1/(4*mm))...
        +(3/4)*(b/a)*(1-1/(4*mm)));
    F_prime(ii_node,ii_node+1)=-(3/8)*(b/a)*((1-1/(4*mm))+eps_ins*(1+1/(4*mm)));
    F_prime(ii_node,ii_node+bound(11)-1)=-(1/8)*(b/a)*(1+1/(4*mm))*eps_ins;
    F_prime(ii_node,ii_node+bound(11))=-((a/b)*(1+1/(2*mm))-(1/4)*(b/a)*(1+1/(4*mm)))*eps_ins;
    F_prime(ii_node,ii_node+bound(11)+1)=-(1/8)*(b/a)*(1+1/(4*mm))*eps_ins;
end

%Specify F_prime at the boundaries.

%The R=0 boundary
j_node = bound(1):bound(2);
F_prime(j_node,:)=0;
for ii_node=bound(1):bound(2)
    F_prime(ii_node,ii_node)=1;
    F_prime(ii_node,ii_node+bound(11))=-1;
end

%The top boundary
j_node = bound(5):bound(8);
F_prime(j_node,:)=0;
for ii_node=bound(5):bound(8)
    if ii_node<bound(6)
        F_prime(ii_node,ii_node)=1;
        F_prime(ii_node,ii_node-bound(11))=-1;
    elseif (ii_node>=bound(6))&(ii_node<=bound(7))
        F_prime(ii_node,ii_node)=1;
    elseif ii_node>bound(7)
        F_prime(ii_node,ii_node)=1;
        F_prime(ii_node,ii_node-bound(11))=-1;
    end
end

%The left and right boundaries

ii_node_l=bound(1);
ii_node_r=bound(2);
for jjj=1:bound(12)
    %Initialization
    F_prime(ii_node_l,:)=rowzero(1,:);
    F_prime(ii_node_r,:)=rowzero(1,:);
    if jjj==1
        F_prime(ii_node_l,ii_node_l)=1;
        F_prime(ii_node_l,ii_node_l+1)=-1/2;
        F_prime(ii_node_l,ii_node_l+bound(11))=-1/2;
        
        F_prime(ii_node_r,ii_node_r)=1;
        F_prime(ii_node_r,ii_node_r-1)=-1/2;
        F_prime(ii_node_r,ii_node_r+bound(11))=-1/2;
        
    elseif (jjj>1)&(jjj<bound(12))
        F_prime(ii_node_l,ii_node_l)=1;
        F_prime(ii_node_l,ii_node_l+1)=-1;
        
        F_prime(ii_node_r,ii_node_r)=1;
        F_prime(ii_node_r,ii_node_r-1)=-1;
        
    elseif jjj==bound(12)
        F_prime(ii_node_l,ii_node_l)=1;
        F_prime(ii_node_l,ii_node_l+1)=-1/2;
        F_prime(ii_node_l,ii_node_l-bound(11))=-1/2;
        
        F_prime(ii_node_r,ii_node_r)=1;
        F_prime(ii_node_r,ii_node_r-1)=-1/2;
        F_prime(ii_node_r,ii_node_r-bound(11))=-1/2;
    end
    ii_node_l=ii_node_l+bound(11);
    ii_node_r=ii_node_r+bound(11);
end
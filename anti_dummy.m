% the anti_dummy function to compute Fn for given Em_cnt and Ne/N1D by the bisectional method.

function [x]=anti_dummy(y,delta,Nband)
Np=length(y);
xup=100*ones(Np,1);
xlow=-100*ones(Np,1);

ii_iter=0; 
while ii_iter<23 
    res=y-dummy((xup+xlow)*0.5,delta,0,Nband);
    error=max(abs(res));
    for ii_np=1:Np
        if res(ii_np)>0
            xlow(ii_np)=(xup(ii_np)+xlow(ii_np))/2;
        else
            xup(ii_np)=(xup(ii_np)+xlow(ii_np))/2;
        end
    end
    ii_iter=ii_iter+1;
end
if (max(abs(res))>1)
    disp('anti_dmmuy exceeds iteration limit');
end
x=(xup+xlow)*0.5;

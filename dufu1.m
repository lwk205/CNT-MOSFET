function y = dufu1(x,zeta,delta)
y = 1./(1+exp(sqrt(x.^2+delta.^2)-zeta)) - 1./(1+exp(sqrt(x.^2+delta.^2)+zeta));

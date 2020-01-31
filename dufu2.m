function y = dufu2(x,zeta,delta)
y = 0.25./cosh((sqrt(x.^2+delta.^2)-zeta)./2).^2 + 0.25./cosh((sqrt(x.^2+delta.^2)+zeta)./2).^2;

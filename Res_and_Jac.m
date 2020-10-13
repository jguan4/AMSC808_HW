function [r,J] = Res_and_Jac(par,xy)
[dim,n] = size(xy); % each point is a column vector
npar = length(par);
[v,W,u] = param(par);
r = zeros(n,1);
J = zeros(n,npar);
for i = 1 : n
    [ri,dri] = res(xy(:,i),v,W,u);
    r(i) = ri;
    J(i,:) = dri';
end

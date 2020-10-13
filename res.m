function [r,dr] = res(x_vec,v,W,u)
% BVP for the Poisson equation is setup here
%
% residual functions r and theor derivatives w.r.t. parameters dr 
% are evaluated in this function
%
% computer diff_operator(Psi(x)) - RHS(x)
% boundary functions
[~,d2A,~,rhs,~] = setup();
% differential operator is d^2/dx^2 + d^2/dy^2
% differential operator applied to A(x,y), the bdry term
d2A_val = d2A(x_vec);
% differential operator applied to B(x,y) = x(1-x)y*phi(x,y,v,W,u)
[phi,phix,phiy,phixx,phiyy,dphi,dphix,dphiy,dphixx,dphiyy] = phi_fun(x_vec,v,W,u);
x = x_vec(1); y = x_vec(2);
d2B = -2*y*phi+2*y*(1-2*x)*phix+2*x*(1-x)*phiy+y*x*(1-x)*(phixx+phiyy);
% residual r = d2A + d2B - RHS
r = d2A_val + d2B - rhs(x_vec(1),x_vec(2));
% derivative of r w.r.t. parameters
dr = -2*y*dphi+2*y*(1-2*x)*dphix+2*x*(1-x)*dphiy+y*x*(1-x)*(dphixx+dphiyy);
end



function [B,d2B,h,rhs,exact_sol] = setup()
% differential operator is d^2/dx^2 + d^2/dy^2
B = @(x,y) 2*y.*sin(pi*x);
% differential operator applied to A(x,y), the bdry term
d2B = @(x) -2*pi^2*x(2)*sin(pi*x(1));
% differential operator applied to B(x,y) = x(1-x)y(NN(x,y,v,W,u)-NN(x,1,v,W,u)-dNNdy(x,1,v,W,u))
h = @(x,y) x.*(1-x).*y;
% right-hand side
rhs = @(x,y)(2-pi^2*y.^2).*sin(pi*x);
exact_sol = @(x,y)y.^2.*sin(pi*x);
end
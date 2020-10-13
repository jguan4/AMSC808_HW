function [phi_val,phix,phiy,phixx,phiyy,dphi,dphix,dphiy,dphixx,dphiyy] = phi_fun(x_vec,v,W,u)
[fun,dfun,d2fun,d3fun,d4fun] = ActivationFun();
x_1 = [x_vec(1);1];
[f,fx,fy,~,fxx,fyy,~,df,dfx,dfy,~,dfxx,dfyy,~] = NN(x_vec,v,W,u,fun,dfun,d2fun,d3fun,d4fun);
[f1,fx1,fy1,fxy1,fxx1,~,fxxy1,df1,dfx1,dfy1,dfxy1,dfxx1,~,dfxxy1] = NN(x_1,v,W,u,fun,dfun,d2fun,d3fun,d4fun);
phi_val = f-f1-fy1;
phix = fx-fx1-fxy1;
phiy = fy;
phixx = fxx-fxx1-fxxy1;
phiyy = fyy;
dphi = df-df1-dfy1;
dphix = dfx-dfx1-dfxy1;
dphiy = dfy;
dphixx = dfxx-dfxx1-dfxxy1;
dphiyy = dfyy;
end
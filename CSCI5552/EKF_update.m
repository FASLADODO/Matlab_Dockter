function[x_hat_plus_k, P_plus_k, res_k, S_k] = EKF_update(x_hat_min_k, P_min_k, z_g_k, sigma_g)
%%Rod Dockter
%%Assignment 2 problem 2
%%function to implement the discrete time kalman filter update equations
%%According to equations 4.30 - 4.37 in section 4.3
%%Assuming
H=1;
z_hat_k=H*x_hat_min_k;
res_k =z_g_k-z_hat_k;
R_k=sigma_g^2; %%R_li = Estimate[sigma*Sigma^T]
S_k=H*P_min_k*H+R_k;
K_k=P_min_k*H*(1/S_k);
x_hat_plus_k=x_hat_min_k+K_k*res_k;
P_plus_k=P_min_k-P_min_k*H*(1/S_k)*H*P_min_k;
end
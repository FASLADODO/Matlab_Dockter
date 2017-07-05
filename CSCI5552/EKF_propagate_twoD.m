function [x_hat_min_kplus,P_min_kplus] = EKF_propagate_twoD(x_hat_plus_k,P_plus_k, theta_hat_est, v_m_k, w_m_k, sigma_v, sigma_w, dt)
%%Problem 3 Assignment 2
%% Rod Dockter
%%equations as follows (assuming constant velocities a=0)
%% According to 4.64 - 4.67
%%x_hat_k+1 = f(x_hat_k, u_k, 0);
%%and
%%P(k+1|k)=phi(k)*P(k|k)*phi(k)^T+G(k)*Q(k)*G(k)^T
%%where now phi_k = df/dx(x_hat, u_k,0)
%%and G_k=df/dw*f(x_hat_k,u_k,0)
x_hat_min_kplus = zeros(3,1);

phi_mat = [1,0,-v_m_k*dt*sin(theta_hat_est);0,1,v_m_k*dt*cos(theta_hat_est);0,0,1];
G_mat=[-dt*cos(theta_hat_est),0;-dt*sin(theta_hat_est),0;0,-dt];

x_hat_min_kplus(3,1) = theta_hat_est + w_m_k*dt;
x_hat_min_kplus(1,1) = x_hat_plus_k(1,1) + v_m_k*dt*cos(theta_hat_est);
x_hat_min_kplus(2,1) = x_hat_plus_k(2,1) + v_m_k*dt*sin(theta_hat_est);

%%note that Q_k=E[W_k*W_k^T|z^k]
%%so from chapter 5 we get

Q_k=[sigma_v^2,sigma_v*sigma_w;sigma_v*sigma_w,sigma_w^2];
P_min_kplus=phi_mat*P_plus_k*transpose(phi_mat)+G_mat*Q_k*transpose(G_mat);

end
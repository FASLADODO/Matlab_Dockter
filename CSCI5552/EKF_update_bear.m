function [x_k_plus, p_k_plus,r_k,s_k] = EKF_update_bear(x_hat_min, p_min, Xldmk, Yldmk, bearing, nl)
%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552

%%EKF update for bearing only
%%This is called by twoD_EKF
%It will return the EKF update for the bearing only exteroceptive
%measurements

%%setting
sigmatheta = 0.01;
z=bearing;
%%getting R_k
R_k=sigmatheta^2*eye(nl);

%%This loop sets up the zhat and H matrix for only the number of landmarks
%%we wish to use as set in twoD_EKF
for i=1:nl
    %%using z=h(x) equation for bearing measurements
    z_hat(i,1) = -x_hat_min(3)+atan2(abs(Yldmk(i)-x_hat_min(2)),abs(Xldmk(i)-x_hat_min(1)));

    %%getting H matrix for bearing
    H_k(i,:) = [(Yldmk(i)-x_hat_min(2))/((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),(-Xldmk(i)+x_hat_min(1))/((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),-1];
end
%%Standard EKF update equations using z and H and returned the updated
%%postion and covariance
r_k = z-z_hat;
s_k = H_k*p_min*H_k'+R_k;
k_k = p_min*H_k'*inv(s_k);
x_k_plus = x_hat_min+k_k*r_k;
p_k_plus=p_min-p_min*H_k'*inv(s_k)*H_k*p_min;

end
function [x_k_plus, p_k_plus,r_k,s_k] = EKF_update_dist(x_hat_min, p_min, Xldmk, Yldmk, distance, nl)
%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552

%%Ekf update for distance measurements only
%%This gets called by the twoD_EKF for the distance only update
%%EKF is done in here for one position step at a time

%%setting
sigmad = 0.01;
z=distance;
%%getting R_k
R_k=sigmad^2*eye(nl);

%%getting z hat and H depending on number of landmarks
%%Only the number of landmarks desired will be used in z and H
for i=1:nl
    %%using z=h(x) equation for distance measurements
    z_hat(i,1) = sqrt((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2);

    %%getting H matrix for distance
    H_k(i,:) = [(-Xldmk(i)+x_hat_min(1))/sqrt((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),(-Yldmk(i)+x_hat_min(2))/sqrt((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),0];
end

%%Using standard EKF equations for the distance measurements and the H
%%matrix as determined.
r_k = z-z_hat;
s_k = H_k*p_min*H_k'+R_k;
k_k = p_min*H_k'*inv(s_k);
x_k_plus = x_hat_min+k_k*r_k;
p_k_plus=p_min-p_min*H_k'*inv(s_k)*H_k*p_min;
%%These x_plus and p_plus values are returned
end
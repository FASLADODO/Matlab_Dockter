function [x_k_plus, p_k_plus,r_k,s_k] = EKF_update_bear_dist(x_hat_min, p_min, Xldmk, Yldmk, bearing, distance, nl)
%%Rod Dockter
%%Homework 3
%%problem 3
%%csci 5552

%%EKF update for bearing and distance
%%This is called by twoD_EKF in order to compute the combined distance and
%%bearing measurements to landmarks. These are computed seperately and the
%%combined into z and H alternatingly. 

%%setting sigmas (both)
sigmad = 0.01;
sigmatheta = 0.01;

%%setting the true measurements into the z matrix (alternating distance and
%%bearings each row)
for j = 1:nl
    z(2*(j-1)+1,1)=distance(j,1);
    z(2*(j-1)+2,1)=bearing(j,1);
end

%%getting R_k (in this case R has to be 6x6 for 3 landmarks and since the
%%noise for bearing and distance are the same I can just use a standard
%%diagonal) this wouldnt work if bearing and distance and different noises
R_k=sigmad^2*eye(2*nl);

%%calculating the zhats and H for the combined distance and bearing, only
%%using the number of landmarks specified in twoD_EKF
for i=1:nl
     %%using z=h(x) equation for distance measurements
    z_hat(2*(i-1)+1,1) = sqrt((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2);

    %%using z=h(x) equation for bearing measurements
    z_hat(2*(i-1)+2,1) = -x_hat_min(3)+atan2(abs(Yldmk(i)-x_hat_min(2)),abs(Xldmk(i)-x_hat_min(1)));

    %%getting H matrix for distance
    H_k(2*(i-1)+1,:) = [(-Xldmk(i)+x_hat_min(1))/sqrt((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),(-Yldmk(i)+x_hat_min(2))/sqrt((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),0];

    %%getting H matrix for bearing
    H_k(2*(i-1)+2,:) = [(Yldmk(i)-x_hat_min(2))/((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),(-Xldmk(i)+x_hat_min(1))/((Xldmk(i)-x_hat_min(1))^2+(Yldmk(i)-x_hat_min(2))^2),-1];
end

%%Again we finally use the standard EKF update equations but now for the
%%larger H, z, res and S matrices
r_k = z-z_hat;
s_k = H_k*p_min*H_k'+R_k;
k_k = p_min*H_k'*inv(s_k);
x_k_plus = x_hat_min+k_k*r_k;
p_k_plus=p_min-p_min*H_k'*inv(s_k)*H_k*p_min;

end
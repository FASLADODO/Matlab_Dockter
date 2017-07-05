function [R1_p_R2_hat, R1_phi_R2_hat, P_xx_inv] = ICP_wls(R1_p_R2_init, R1_phi_R2_init, scan1, scan2, sigma_r, sigma_th)
%%Rod Dockter
%%Assignment 2, Problem 4
%%This program implemenets the weighted least squres estimation method,
%%operating on the measured r and phi from the rws scans of the two
%%objects.
%%To do this I pass in the values of r and phi from the scans as well as
%%the matched indices from the ICP_match function in order to know the non-zero
%%elements of the scans.
%%Then using the r and phi values, convert them to cartesian
%%coordinates at all the non-zero matches. Then using these new x and y
%%coordinates from each scan, build our xhat (3x1) and z (2x1) vector.
%%Then I build the summation matrices according to chapter 3, equations
%%3.43-45. As well as using the equations for problem 5 of homework 1 (pg.8)
%%which are the weighted least squares estimate for the localization
%%problem:
%%xhatk+1 = xhatk + sum(Hi'*Ri^-1*Hi)^-1*sum(Hi;*Ri^-1*(R1_zi+R1_p_R2+R1_C(phi)_R2*R2_z_i))
%%For ease of coding I will call:
%%matrix1 = sum(Hi'*Ri^1*Hi) since it will get used a lot, and
%%matrix2 = sum(Hi;*Ri^-1*(R1_zi+R1_p_R2+R1_C(phi)_R2*R2_z_i))
%%Where Hi = [I J*R1_C(phi)_R2*R2_zi] and Ri = R1_R_Li + R1_C_R2*R2_R_Li*R1_C_R2
%%Where j = [0,-1;1,0] (56) and z = [x;y] and R1_p_R2 = [x;y;phi] 
%%And lastly Ri comes from the get Ri function
%%This yields our value for the new x hat propogation
%%Then we get covariances from equation (3.46) in chapter 3 as 
%%pxx = sum(Hi'*Ri^-1*Hi)^-1
%%Forgive my mathematical innacuracies in writing these equations out, they
%%are all in the notes anyway.

%%Getting matches
[match1, match2] = ICP_match(R1_p_R2_init,R1_phi_R2_init,scan1,scan2);

%%Initializing extra matrices for summation and the j and C matrices
matrix1 = zeros(3,3);
matrix2 = zeros(3,1);
j = [0 -1; 1 0];
C_phi = [cos(R1_phi_R2_init),-sin(R1_phi_R2_init);sin(R1_phi_R2_init),cos(R1_phi_R2_init)];
%%initializing the xhat
xhat_vector = [R1_p_R2_init; R1_phi_R2_init];

%%Looping through the summation of the matrices
for q=1:length(match1)
    %%Calling the get_ri function
    R1_Ri = get_Ri(scan1.rm(match1(q)), scan1.theta(match1(q)), sigma_r, sigma_th);
    R2_Ri = get_Ri(scan2.rm(match2(q)), scan2.theta(match2(q)), sigma_r, sigma_th);
    %%Computing the overall Ri matrix
    R_i = R1_Ri +C_phi*R2_Ri*C_phi';
    %%Converting to cartesian coordinates
    [x_scan1, y_scan1] = pol2cart(scan1.theta(match1(q)), scan1.rm(match1(q)));
    [x_scan2, y_scan2] = pol2cart(scan2.theta(match2(q)), scan2.rm(match2(q)));
    %%Building the z vector
    scan1_vector = [x_scan1; y_scan1];
    scan2_vector = [x_scan2; y_scan2];
    %%building the H matrix
    H_i = [eye(2),j*C_phi*scan2_vector];
    %%Updating the summation of the two matrices
    %%Unfortunately matlab dislikes Ri\H for different sizes, forcing the
    %%use of inv()
    matrix1 = matrix1 + H_i'*inv(R_i)*H_i;
    matrix2 = matrix2 + H_i'*inv(R_i)*(scan1_vector - R1_p_R2_init - C_phi*scan2_vector);
end
%%WLS for the location method using the two defined sub matrices
xhat_vector = xhat_vector + inv(matrix1)*matrix2;

R1_p_R2_hat = xhat_vector(1:2,1);
R1_phi_R2_hat = xhat_vector(3,1);
%%Calculate the non-inverted component of covariance according to equation
%%3.46
P_xx_inv = matrix1;

end
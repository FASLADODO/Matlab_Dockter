function [x_hat_p, P_p, res_update, S_update, Nc]= SLAM_EKF_Update(x_hat_m, P_m, z, Nc, R)
%%Rod Dockter
%%EKF update equations
%%CSCI 5552
%%ported from c++ code for the project

% z_hat(k+1|k) = h[x_hat(k+1|k),u(k),0]
% r(k+1|k) = z(k+1) - z_hat(k1|k)
% S(k+1|k) = H(k+1)*P(k+1|k)*H(k+1)' + R(k+1)
% K(k+1|k) = P(k+1|k)*H(k+1)'*inv[S(k+1|k)]
% x_hat(k+1|k+1) = x_hat(k+1|k) + K(k+1|k)*r(k+1|k)
% P(k+1|k+1) = P(k+1|k) - P(k+1|k)*H(k+1)'*inv[S(k+1|k)]*H(k+1)*P(k+1|k)

%%getting x - y coordinates of landmarks from d and theta
for i1 = 1:(length(z)/2)
    zpnt((2*i1-1):(2*i1),1) = [z(2*i1-1)*cos(z(2*i1));
                             z(2*i1-1)*sin(z(2*i1))];
end

%Jacobian
J = [0 -1; 1 0];

%rotaitonal matrix
C = [cos(x_hat_m(3,1)),-sin(x_hat_m(3,1));...
     sin(x_hat_m(3,1)),cos(x_hat_m(3,1))];

 %%setting upper and lower bounds for updating and intiliazing respectively
 %%with mahalanobis distance
gammaUpper = 50;
gammaLower = 6;

%%Initializing update and initialization values
nUpdate = 0;
nInit = 0;

% Iterating through measurements
for i1 = 1:(length(zpnt)/2)
    if (Nc > 0)
        % looking at all previous corners
        for i2 = 1:Nc
            H = zeros(2,3+2*Nc);
            H(1:2,1:2) = -C';
            H(1:2,3) = -C'*J*(x_hat_m((2*i2+2):(2*i2+3),1)-x_hat_m(1:2,1));
            H(1:2,(2*i2+2):(2*i2+3)) = C';
            
            %getting z estimate and residual
            z_hat = C'*(x_hat_m((2*i2+2):(2*i2+3),1)-x_hat_m(1:2,1));
            r = zpnt((2*i1-1):(2*i1),1) - z_hat;
            
            G = [cos(z(2*i1)), -z(2*i1-1)*sin(z(2*i1));...
                 sin(z(2*i1)), z(2*i1-1)*cos(z(2*i1))];
            
            % Compute the covariance <S> of the residual
            S = H*P_m*H'+G*R*G';
            
            % Compute the mahalanobis distance <gamma>
            gamma(i2) = r'*inv(S)*r;
        end
        % calculating a minimum mahalanobis distance
        [gammaMin gammaID] = min(gamma);
    end
    
    % Verify if corner is from landmark in state
    if ((Nc > 0) && (gammaMin < gammaLower))
        Lm_ID(nUpdate+1) = gammaID;
        zpnt_update((2*nUpdate+1):(2*nUpdate+2),1) = zpnt((2*i1-1):(2*i1),1);
        z_update((2*nUpdate+1):(2*nUpdate+2),1) = z((2*i1-1):(2*i1),1);
        nUpdate = nUpdate+1;
    elseif (( Nc == 0) || (gammaMin > gammaUpper))
        zpnt_init((2*nInit+1):(2*nInit+2),1) = zpnt((2*i1-1):(2*i1),1);
        z_init((2*nInit+1):(2*nInit+2),1) = z((2*i1-1):(2*i1),1);
        nInit = nInit + 1;
    end
end

% If previous landmark was found again perform an update
if (nUpdate > 0)
    H_update = zeros(2*nUpdate,3+2*Nc);
    res_update = zeros(2*nUpdate,1);
    R_update = zeros(2*nUpdate,2*nUpdate);
        
    % Building EKF H and R and S for corners
    for k = 1:nUpdate
        H_update((2*k-1):(2*k),1:2) = -C';
        H_update((2*k-1):(2*k),3) = -C'*J*(x_hat_m((2*Lm_ID(k)+2):(2*Lm_ID(k)+3),1)-x_hat_m(1:2,1));
        H_update((2*k-1):(2*k),(2*Lm_ID(k)+2):(2*Lm_ID(k)+3)) = C';
        
        % Compute the measurement estimate <z_hat> and corresponding residual <r>
        z_hat = C'*(x_hat_m((2*Lm_ID(k)+2):(2*Lm_ID(k)+3),1)-x_hat_m(1:2,1));
        
        res_update((2*k-1):(2*k),1) = zpnt_update((2*k-1):(2*k),1) - z_hat;
        
        G = [cos(z_update(2*k)), -z_update(2*k-1)*sin(z_update(2*k));...
        	sin(z_update(2*k)), z_update(2*k-1)*cos(z_update(2*k))];
        
        % Add the variance of the measurements to the measurement covariance matrix <R>
        R_update((2*k-1):(2*k),(2*k-1):(2*k)) = G*R*G';
    end
    
    % calculate covariance of the residual
    S_update = H_update*P_m*H_update' + R_update;
    
    % calculate the kalman gain
    K_update = P_m*H_update'*inv(S_update);
    
    % Updating the state
    x_hat_p = x_hat_m + K_update*res_update;
    
    % Updating state covariance 
    P_p = (eye(3+2*Nc)-K_update*H_update)*P_m*(eye(3+2*Nc)-K_update*H_update)' + K_update*R_update*K_update';
    
end

% If new landmarks were identified, initilialize
if (nInit > 0)
    
    Hr = zeros(2,3);
            
    if (nUpdate > 0)
        C = [cos(x_hat_p(3,1)),-sin(x_hat_p(3,1));...
             sin(x_hat_p(3,1)),cos(x_hat_p(3,1))];
        
        Hc = C';
        Hr(1:2,1:2) = -C';
        
        % Iterating new landmark measurements
        for k = 1:nInit
            % Add landmarks to state using update
            x_hat_p((2*(Nc+k)+2):(2*(Nc+k)+3),1) = x_hat_p(1:2,1) + C*zpnt_init((2*k-1):(2*k),1);
                        
            Hr(1:2,3) = -C'*J*(x_hat_p((2*(Nc+k-1)+4):(2*(Nc+k-1)+5),1)-x_hat_p(1:2,1));
            
            G = [cos(z_init(2*k)), -z_init(2*k-1)*sin(z_init(2*k));...
                sin(z_init(2*k)), z_init(2*k-1)*cos(z_init(2*k))];
        
            P_p = [P_p,-P_p(1:(2*(Nc+k)+1),1:3)*Hr'*Hc;...
                    -Hc'*Hr*P_p(1:3,1:(2*(Nc+k)+1)),Hc'*(Hr*P_p(1:3,1:3)*Hr' + G*R*G')*Hc];
        end
    else
        C = [cos(x_hat_m(3,1)),-sin(x_hat_m(3,1));...
             sin(x_hat_m(3,1)),cos(x_hat_m(3,1))];
        
        Hc = C';
        Hr(1:2,1:2) = -C';
        
        x_hat_p = x_hat_m;
        P_p = P_m;

        % Iterate through all new landmark measurements
        for k = 1:nInit
            % Adding landmarks to state using update estimate
            x_hat_p((2*(Nc+k)+2):(2*(Nc+k)+3),1) = x_hat_p(1:2,1) + C*zpnt_init((2*k-1):(2*k),1);
                        
            Hr(1:2,3) = -C'*J*(x_hat_p((2*(Nc+k-1)+4):(2*(Nc+k-1)+5),1)-x_hat_p(1:2,1));
            
            G = [cos(z_init(2*k)), -z_init(2*k-1)*sin(z_init(2*k));...
                sin(z_init(2*k)), z_init(2*k-1)*cos(z_init(2*k))];
        
            P_p = [P_p,-P_p(1:(2*(Nc+k)+1),1:3)*Hr'*Hc;...
                    -Hc'*Hr*P_p(1:3,1:(2*(Nc+k)+1)),Hc'*(Hr*P_p(1:3,1:3)*Hr' + G*R*G')*Hc];
        end
    end
    %%increment number of landmarks
    Nc = Nc + nInit;
end

%Only resturning residual if we got an update
if (nUpdate == 0)
   res_update = [];
   S_update = [];
end

%%if no update, propogate state
if ((nUpdate == 0) && (nInit ==0))
    x_hat_p = x_hat_m;
    P_p = P_m;
end
end
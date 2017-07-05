%make a real slam ekf
%http://www.iri.upc.edu/people/jsola/JoanSola/objectes/curs_SLAM/SLAM2D/SLAM%20course.pdf

clear all

N = 20; %number of landmarks

%set initial landmarks
W = [randn(1,N)*5 + 5; randn(1,N)*5 + 5 ];


%set initial robot locations
% R: robot pose [x ; y ; alpha]
R = [0;-2;0.4];
goalloc = [10;10;0];

%Define motion parameters
deltaT = 0.1; %every 100 ms
velocity = 1.0; %m/s
U_0 = [velocity; 0.1];

%landmark measurements
Y = zeros(2,N);

% I.2 ESTIMATOR
% Map: Gaussian {x,P}
% x: state vector's mean
r = [1:3]; %robot state indices
x(r,:) = R;
% P: state vector's covariance
P = 0.01*randn(numel(x),numel(x));
% Jacobians (2 x (3+N) ) it will grow
H = zeros(2,numel(x));

% System noise: Gaussian {0,Q}
q = [0.001;0.00002]; % amplitude or standard deviation
Q = diag(q.^2); % covariances matrix

% Measurement noise: Gaussian {0,S}
s = [0.01;0.1*pi/180]; % amplitude or standard deviation
S = diag(s.^2); % covariances matrix


%landmark indices
knownlandmarks = [];
state = 0; % 1 for update, 2 for add to state

%initial setup
figure(10);
h1 = scatter(W(1,:), W(2,:), 30, 'r*');
hold on
%actual location
h2 = quiver(x(1),x(2),cos(x(3)),sin(x(3)),'g','LineWidth',3);
hold off
legend([h1(1),h2(1)],'Landmarks','robot pose');


%%

ff = figure(2);
axis([-20 20 -20 20]) % set axes limits
axis square % set 1:1 aspect ratio

updatecnt = 0;
updatesteps = 1; %update every two propogates
idxli = 1; %initialize

fprintf('running simulation ...');
while(1)
    %stopping conditions 
    if(norm(goalloc(1:2) - R(1:2)) < 2)
       break; 
    end
    
    if(norm(R(1:2) - x(1:2)) > 1)
       shit = 1; 
    end
    
    %%%%%%%%%%%%%%%%%%%% FOR SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %see if we gotta turn
    deltatheta = atan2(goalloc(2) - R(2), goalloc(1) - R(1)) - R(3);
    U = [deltaT; deltatheta] .* U_0;
    
    %update true position
    n = q .* randn(2,1); % perturbation vector
    [R] = odometryModel(R, U, zeros(2,1) ); % true position

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%SLAM EKF START HERE%%%%%%%%%%%%%%%%%%%%%%
    % b. Prediction ?? robot motion and jacobian of prediction model
    [x(r), A.r, A.n] = odometryModel(x(r), U, n); % Estimator perturbed with n\
    
    %propogate covariance
    P = propogateCovariance(P,A,Q);
    
    % All observations
    for i = 1:N % i: landmark index
        v = s .* randn(2,1); % measurement noise
        %Y = [range;bearing] in true robot frame
        Y(:,i) = observe(R, W(:,i)) + v;
    end


    if(updatecnt >= updatesteps)
        %only update if we have done enough time steps
        updatecnt = 0;
    
        %%%%%%%%%%%%%%%%%%%%%SIMULATE LANDMARKS %%%%%%%%%%%%%%%%%%%%%%%%%%
        idxli = getNewLandmark(Y);
        
        knownlandmarks
        %%%%%%%%%%%%%%%%%%%%%CHECK IF KNOWN LANDMARKS %%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_known = 1;
        %check if we've seen this one before (ignore malhalonobisfor now)
        if(any(knownlandmarks==idxli) )
            state = 1;
            idx_known = find(knownlandmarks == idxli);
        else
            knownlandmarks = [knownlandmarks; idxli];
            state = 2;
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%SLAM EKF UPDATE%%%%%%%%%%%%%%%%%%%%%%

        %Now get noisy measurement signals (this would normally just come form the
        %laser scan
        if(state == 1)%we're seen this before
            %prior landmark location
            idx_state = 2 + idx_known*2;
            l = [idx_state:idx_state+1]; %landmark state indices
            
            %measurement model and jacobians h=[range; bearing] in current
            %robot frame
            [h, E.r, E.l] = observe(x(r), x(l) ); 
            
            %inovation gaussian
            Yi = Y(:,idxli);
            z = (Yi - h) %residual
            [z(2)] = boundAngle(z(2));

            %jacobian update from landmark
            H = updateJacobian(H,E,idx_state);

            %Kalman gains
            [K, gotime] = computeKalmanGain(P,H,S);
            if(gotime == 1)
                %update state
                x = updateStateVector(x,K,z);
                %update covariance
                [P] = updateCovariance(P,K,H,S);
            end
        elseif(state == 2)%i've no memory of this place
            %measurement
            Yi = Y(:,idxli);

            % initialization
            [newScan, L.r, L.y] = invObserve(x(r), Yi);
            
            %Add to the state vector
            x = [x ; newScan];

            %Add new landmarks
            P = addNewLandmarkCovariance(P,L,S);
            %add to jacobian
            H = [H, zeros(2,2)];
        end
        
    end

    clf(ff)

    %landmarks
    h1 = scatter(W(1,:), W(2,:), 30, 'r*');
    hold on
    h2 = scatter(W(1,idxli), W(2,idxli), 30, 'b*');
    hold on
    h3 = quiver(R(1),R(2), Y(1,idxli)*cos(Y(2,idxli)+R(3)), Y(1,idxli)*sin(Y(2,idxli)+R(3)),'c');
    hold on

    %actual location
    h4 = quiver(R(1),R(2),cos(R(3)),sin(R(3)),'g','LineWidth',3);
    %odom location
    h5 = quiver(x(1),x(2),cos(x(3)),sin(x(3)),'m','LineWidth',3);
    hold on
    r2 = [1:2];
    [xx,yy] = cov2elli(x(r2),P(r2,r2),3,16); % x? and y? coordinates of contour
    plot(xx,yy);
    hold off
    legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'True Landmark','Seen landmarks','Measurements','actual Location','est location');

    pause(0.1);
   
    
    %increment counter
    updatecnt = updatecnt + 1;
end
fprintf('done\n');
knownlandmarks




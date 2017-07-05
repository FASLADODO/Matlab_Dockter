%Example code using robot movement and landmark measurements to localize
%robot

%% Create some landmarks and intialize robot location

NL = 10;
landmarks = [ randn(NL,1)*3 + 3,randn(NL,1)*3 + 3] ;


location_init = randn(1,2);
location_goal = randn(1,2) + [3,3];
location_true = location_init;

% are measurements noise
dist_noise = 0.2;
bear_noise = 0.005;


%% Now create particles and bootstrap

NP = 50;
SigmaP = 2; %particle covariance
bw = 1; %bandwidth
epochs = 50; %how many times to run it
scale = 100;
res_thresh = 0.8;
shift = 0.5;

%create some damn particles
particles_init = randn(NP,2)*SigmaP + repmat(location_true,NP,1);


figure
%loop through however many timesteps we're doing
for tt = 1:epochs
    %what are the current particles
    particles = particles_init;
    
    %update location
    stepmove = (location_goal - location_init);
    location_true = location_init + (stepmove * (tt/epochs) );
    rloc = location_true; %robot location
    
    %get measurements to landmarks from robot location
    r_measure = []; %for measurements
    for ll = 1:NL
        lloc = landmarks(ll,:);
        %get distance and angle to it
        step = lloc - rloc;
        dist = norm(step) + randn(1)*dist_noise;
        angle = atan2(step(2),step(1)) + randn(1)*bear_noise;
        r_measure(ll,:) = [dist,angle];
    end
    
    %Get measurements from all particles to all landmarks
    particle_weights = [];
    for pp = 1:NP
        ploc = particles(pp,:); %particle location
        p_measure = [];
        for ll = 1:NL
            lloc = landmarks(ll,:);
            %get distance and angle to it
            step = lloc - ploc;
            dist = norm(step) + randn(1)*dist_noise;
            angle = atan2(step(2),step(1)) + randn(1)*bear_noise;
            p_measure(ll,:) = [dist,angle];
        end
        
        %compute particle weight
        error = sum(mean(abs(r_measure-p_measure)));
        w = exp(-(bw*error)^2); %radial basis error
        particle_weights(pp,:) = w; %save particle weights
    end
    %make weights sum to one
    particle_weights = particle_weights./sum(particle_weights);
    
    %RESAMPLE STEP (THERE ARE LOTS OF WAY TO DO THIS)
    %get the best estimate
    [best_weight,idn] = max(particle_weights);
    best_estimate = particles(idn,:);
    %get ratio from best and decide if we resample
    ratio_import = particle_weights ./ best_weight;
    keep_particles = ratio_import > res_thresh;
    
    %decide on new variance
    sigma_res = (1 - (1/best_weight)) + shift;
    particles_new = randn(NP,2)*sigma_res + repmat(best_estimate,NP,1);
    for pp = 1:NP
       %decide if we keep or not
       if(keep_particles(pp))
            particles_init(pp,:) = particles(pp,:);
       else
            particles_init(pp,:) = particles_new(pp,:);
       end
    end
    
    
    % HASHTAG PLOTZ
    clf
    scatter(landmarks(:,1),landmarks(:,2),'kx')
    hold on
    scatter(particles(:,1),particles(:,2),particle_weights*scale,'bo')
    hold on
    scatter(location_true(:,1),location_true(:,2),'r+')
    hold on
    scatter(location_goal(:,1),location_goal(:,2),'c+')
    hold on
    scatter(best_estimate(:,1),best_estimate(:,2),'go')
    hold off
    legend('landmarks','particles','robot','goal','best est')
    pause(0.5);
    
end




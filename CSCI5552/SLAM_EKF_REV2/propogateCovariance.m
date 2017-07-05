function P = propogateCovariance(P,A,Q)
%update  the robot covariance and existing landmark covariances
%P = existing covariance 3 + 2*n
%A = prediction jacobian
%Q = process noise
%Page 37
%https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-412j-cognitive-robotics-spring-2005/projects/1aslam_blas_repo.pdf


    %number of landmarks
    n = size(P,1);
    nli = (n - 3) / 2;
    r = [1:3]; %robot indices
    
    %Covariance: robot <-> robot
    PR = A.r*P(r,r)*A.r + A.n*Q*A.n';
    P(r,r) = PR;
    
    %update each of the landmarks 3+2*n jacobians as well
    for ii = 1:nli
       %existing landmark covariances
       idx = 2 + 2*ii;
       m = [idx:idx+1];
       
       %Covariance: landmarks <-> robot
       Pri = A.r*P(r,m);
       P(r,m) = Pri;
       P(m,r) = Pri';
    end
end
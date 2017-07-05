function H = updateJacobian(H,E,idxekf)
    %update H
    %E : E.r and E.l (jacobians from observation)
    %idxekf index of landmark in state matrix
    %Page 39
    %https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-412j-cognitive-robotics-spring-2005/projects/1aslam_blas_repo.pdf
   
    %indices 
    all = [1:2];
    r = [1:3];
    l = [idxekf:idxekf+1];
    
    %first 3 columns are just robot state
    H(all,r) = E.r;
    
    %we also stash -A, -B, -D, -E in the corresponding columns for the
    %current landmark
    H(all,l) = E.l;
end
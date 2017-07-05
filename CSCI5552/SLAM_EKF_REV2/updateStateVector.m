function X = updateStateVector(X,K,res)
    %X=state vector
    %K kalman gain
    %res: residual error (z-h)
    %Page 40
    %https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-412j-cognitive-robotics-spring-2005/projects/1aslam_blas_repo.pdf

    X = X + K*res;
end
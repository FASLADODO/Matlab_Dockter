function [radius, center] = bounding_sphere(data,alpha_val)
%bounding_sphere(data,alpha_val); Find a bounding sphere for 3D data given a confidence.
%ie: alpha_val = 0.95 would return a sphere containing 95% of data
%Returns the radius and center

    % Calculate the eigenvectors and eigenvalues
    covariance = cov(data);
    [eigenvec, eigenval ] = eig(covariance);

    % Get the mean eig
    mean_eigenval = mean(diag(eigenval));

    % Get the coordinates of the data mean
    center = mean(data);

    % Get the 95% confidence interval error ellipse
    chisquare_val=sqrt(chi2inv(alpha_val, 3));

    %get the radius
    radius=chisquare_val*sqrt(mean_eigenval);

end
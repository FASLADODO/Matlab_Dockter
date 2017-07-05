function h = MinBoundEllipse_Plot(A,C,n,plotarg)
%Plot ellipse in 2D only
%A: covariance matrix from MivVolEllipse()
%C: centroid from  MivVolEllipse()
%n: number of points on ellipse to plot
%returns h a handle to plot;

if(nargin == 3)
    plotarg = 'c-';
end

% Get eigenvalues and rotation matrix
[~, Q, V] = svd(A);

%get semimajor axis radius
r = 1./sqrt(diag(Q));
rx = r(1);
ry = r(2);

%interpolate some angles
theta = linspace(0,2*pi,n)';

%get regular ellipse coordinates
ellps = [rx*cos(theta), ry*sin(theta)];

%rotate ellipse by rotation matrix V
ellps = (( V * ellps') + repmat(C,1,n))';

%plot and save handle
h = plot(ellps(:,1),ellps(:,2),plotarg);

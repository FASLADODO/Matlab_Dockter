function h = plot_error_ellipse(x_hat,P,varargin)
npts = 100;
th = linspace(0, 2*pi, npts);

% get radius, two different possible methods...
%r = sqrt(chi2inv(0.9973,2)); % here you can specify the percentile
r = 3; % here you can specify the sigma

[Evec, Eval] = eig(P);

pts = repmat(x_hat,1,npts) + r*Evec*sqrt(Eval)*[cos(th);sin(th)];

h1 = plot(pts(1,:), pts(2,:),varargin{:});

if nargout
    h=h1;
end
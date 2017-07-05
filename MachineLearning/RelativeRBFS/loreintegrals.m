% http://blogs.mathworks.com/loren/2006/04/26/two-dimensional-integration-over-a-general-domain/

a = 4;
b = 3;

func = @(x,y) ones(size(x)).*y;

dblquad(func, 0, a, 0, b)


x = 0:a;
y = 0:b;

[XX,YY] = meshgrid(x,y);

data = [XX(:),YY(:)];

z = func(data(:,1),data(:,2));

figure
scatter3(data(:,1),data(:,2),z)

scale = ScaleRBF(data,z)


%% compare integral2 to rods fancy estimate
fun = @(x,y) sin(x) +cos(y)
funl = @(x) sin(x(:,1)) +cos(x(:,2))

limits = [0,0;pi/2,pi/2]

q = integral2(fun,0,pi/2,0,pi/2)

integral_result = IntegralND(funl,limits,[100,100])




%% more complicated

func = @(x,y) x.*y;

a = 4;
b = 3;

tol = [];
qnumeric = dblquad(func, 0, a, 0, b, tol, @quadl)
qexact = (a^2 * b^2) / 4


%rods estimate
x = 0:0.1:a;
y = 0:0.1:b;

[XX,YY] = meshgrid(x,y);

data = [XX(:),YY(:)];

z = func(data(:,1),data(:,2));

figure
scatter3(data(:,1),data(:,2),z)

qestimate = ScaleRBF(data,z)


%% gaussian
func = @(x,y) exp(-(x.^2+y.^2) );


tol = [];
qnumeric = dblquad(func, -4, 4, -3, 3, tol, @quadl)

%rods estimate
x = -4:0.1:4;
y = -3:0.1:3;

[XX,YY] = meshgrid(x,y);

data = [XX(:),YY(:)];

z = func(data(:,1),data(:,2));

figure
scatter3(data(:,1),data(:,2),z)

qestimate = ScaleRBF(data,z)

newintegral = dblquad(func, -4, 4, -3, 3, tol, @quadl)/qestimate


%% Try

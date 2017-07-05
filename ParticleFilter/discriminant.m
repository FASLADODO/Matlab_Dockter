%http://www.mathworks.com/help/releases/R2013b/stats/discriminant-analysis.html#bs2oujj

%Load the data into meas and species

load fisheriris;

PL = meas(:,3);
PW = meas(:,4);

h1 = gscatter(PL,PW,species,'krb','ov^',[],'off');
legend('Setosa','Versicolor','Virginica','Location','best')
hold on

%Create a linear discriminant analysis classifier:
X = [PL,PW];
linclass = ClassificationDiscriminant.fit(X,species);

%Create a quadratic classifier:
cqs = ClassificationDiscriminant.fit(X,species,...
    'DiscrimType','quadratic');

K = linclass.Coeffs(2,3).Const; % First retrieve the coefficients for the linear
L = linclass.Coeffs(2,3).Linear;% boundary between the second and third classes
                           % (versicolor and virginica).

% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2 = ezplot(f,[.9 7.1 0 2.5]);
set(h2,'Color','r','LineWidth',2)

hold on

% Now, retrieve the coefficients for the linear boundary between the first
% and second classes (setosa and versicolor).
K = linclass.Coeffs(1,2).Const;
L = linclass.Coeffs(1,2).Linear;

% Plot the curve K + [x1,x2]*L  = 0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h3 = ezplot(f,[.9 7.1 0 2.5]);
set(h3,'Color','k','LineWidth',2)
axis([.9 7.1 0 2.5])
xlabel('Petal Length')
ylabel('Petal Width')
title('{\bf Linear Classification with Fisher Training Data}')

%% Quadtratic

cqs = ClassificationDiscriminant.fit(X,species,...
    'DiscrimType','quadratic');
%Plot the classification boundaries similarly.

%delete(h2); delete(h3) % First, remove the linear boundaries from the plot.

% Now, retrieve the coefficients for the quadratic boundary between the
% second and third classes (versicolor and virginica).
K = cqs.Coeffs(2,3).Const;
L = cqs.Coeffs(2,3).Linear;
Q = cqs.Coeffs(2,3).Quadratic;

% Plot the curve K + [x1,x2]*L + [x1,x2]*Q*[x1,x2]' = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;
h2 = ezplot(f,[.9 7.1 0 2.5]);
set(h2,'Color','r','LineWidth',2)

% Now, retrieve the coefficients for the quadratic boundary between the
% first and second classes (setosa and versicolor).
K = cqs.Coeffs(1,2).Const;
L = cqs.Coeffs(1,2).Linear;
Q = cqs.Coeffs(1,2).Quadratic;

% Plot the curve K + [x1,x2]*L + [x1,x2]*Q*[x1,x2]'=0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;
h3 = ezplot(f,[.9 7.1 0 1.02]); % Plot the relevant portion of the curve.

set(h3,'Color','k','LineWidth',2)
axis([.9 7.1 0 2.5])
xlabel('Petal Length')
ylabel('Petal Width')
title('{\bf Quadratic Classification with Fisher Training Data}')
hold off



%% 2014 version

%Create a linear discriminant analysis classifier:
linclass = fitcdiscr(meas,species);

%Classify an iris with average measurements:
meanmeas = mean(meas);
meanclass = predict(linclass,meanmeas)

meanclass = 'versicolor'

%Create a quadratic classifier:
quadclass = fitcdiscr(meas,species,'discrimType','quadratic');

%Classify an iris with average measurements using the quadratic classifier:
meanclassq = predict(quadclass,meanmeas)

meanclassq = 'versicolor'
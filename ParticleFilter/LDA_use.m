X = [randn(10,2); randn(15,2) + 1; randn(10,2) + 4];
Y = [zeros(10,1); ones(15,1); ones(10,1)*2];

figure(1)
scatter(X(1:10,1),X(1:10,2),'r')
hold on
scatter(X(11:25,1),X(11:25,2),'g')
hold on
scatter(X(26:35,1),X(26:35,2),'c')
hold off
title('actual classes')

% Calculate linear discriminant coefficients
W = LDA(X,Y);

% Calulcate linear scores for training data
L = [ones(25,1), X] * W';

% Calculate class probabilities
P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);

figure(2)
for qq = 1:length(P)
    if(P(qq,1) > P(qq,2))
        scatter(X(qq,1),X(qq,2),'r')
    else
        scatter(X(qq,1),X(qq,2),'g')
    end
    hold on
end

hold off
title('predicted classes')

%% LDA example

clear all

load fisheriris

PL = meas(:,3);
PW = meas(:,4);

scatter(PL,PW)

%%

h1 = gscatter(PL,PW,species,'krb','ov^',[],'off');
set(h1, 'LineWidth', 2);
%legend('Setosa','Versicolor','Virginica','Location','best')
hold on

%creating linear discriminator
X = [PL,PW];
cls = fitcdiscr(X,species);

K = cls.Coeffs(2,3).Const; % First retrieve the coefficients for the linear
L = cls.Coeffs(2,3).Linear;% boundary between the second and third classes
                           % (versicolor and virginica).

% Plot the curve K + [x,y]*L  = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h2 = ezplot(f,[.9 7.1 0 2.5]);
set(h2, 'Color', 'r', 'LineWidth', 2);

% Now, retrieve the coefficients for the linear boundary between the first
% and second classes (setosa and versicolor).
K = cls.Coeffs(1,2).Const;
L = cls.Coeffs(1,2).Linear;

% Plot the curve K + [x1,x2]*L  = 0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
h3 = ezplot(f,[.9 7.1 0 2.5]);
set(h3, 'Color', 'k', 'LineWidth', 2);
axis([.9 7.1 0 2.5])
xlabel('Petal Length')
ylabel('Petal Width')
title('{ Linear Classification with Fisher Training Data}')

%% QDA example

clear all

load fisheriris

PL = meas(:,1);
PW = meas(:,2);

h1 = gscatter(PL,PW,species,'krb','ov^',[],'off');
set(h1, 'LineWidth', 2);
%legend('Setosa','Versicolor','Virginica','Location','best')
hold on

%creating linear discriminator
X = [PL,PW];

cqs = fitcdiscr(X,species,...
    'DiscrimType','quadratic');


% Now, retrieve the coefficients for the quadratic boundary between the
% second and third classes (versicolor and virginica).
K = cqs.Coeffs(2,3).Const;
L = cqs.Coeffs(2,3).Linear;
Q = cqs.Coeffs(2,3).Quadratic;

% Plot the curve K + [x1,x2]*L + [x1,x2]*Q*[x1,x2]' = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;
h2 = ezplot(f,[.9 7.1 0 2.5]);
set(h2, 'Color', 'r', 'LineWidth', 2);

% Now, retrieve the coefficients for the quadratic boundary between the
% first and second classes (setosa and versicolor).
K = cqs.Coeffs(1,2).Const;
L = cqs.Coeffs(1,2).Linear;
Q = cqs.Coeffs(1,2).Quadratic;

% Plot the curve K + [x1,x2]*L + [x1,x2]*Q*[x1,x2]'=0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;
h3 = ezplot(f,[.9 7.1 0 1.02]); % Plot the relevant portion of the curve.
set(h3, 'Color', 'k', 'LineWidth', 2);
axis([.9 7.1 0 2.5])
xlabel('Petal Length')
ylabel('Petal Width')
title('{ Quadratic Classification with Fisher Training Data}')
hold off
                           

%%
n = 1000;
p = [(n-3)/n,1/n,1/n,1/n]
e = EntropyCalc(p)
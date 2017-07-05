% Homework 3, Problem 2
% Analytical Dyad program for 3PP motion generation
% Rod Dockter, Oct. 2014

% setup
clear all
deg2rad = (pi/180);

%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%
% Specify Beta (angle between input positions) (deg)
beta2 = 71;
beta3 = 150;

% Precision positions real+i*im (x,y)
PP1 = 0 + 0*i;
PP2 = -200 + 600*i;
PP3 = -1000 + 800*i;

% Precision angles (deg)
theta1 = 0;
theta2 = 40;
theta3 = 90;



%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%
% convert to radians because matlab
theta1 = theta1 * deg2rad;
theta2 = theta2 * deg2rad;
theta3 = theta3 * deg2rad;
beta2 = beta2 * deg2rad;
beta3 = beta3 * deg2rad;

% position differences (r*e^(i*theta) form )
delta2 = abs(PP2 - PP1)*exp(i*angle(PP2 - PP1));
delta3 = abs(PP3 - PP1)*exp(i*angle(PP3 - PP1));

%angle differences
alpha2 = theta2-theta1;
alpha3 = theta3-theta1;


% Cramers Rule Matrices
denominator = [ exp(i*beta2) - 1, exp(i*alpha2) - 1;
                exp(i*beta3) - 1, exp(i*alpha3) - 1 ];
numerator1 = [ delta2, exp(i*alpha2) - 1;
                delta3, exp(i*alpha3) - 1 ];
numerator2 = [ exp(i*beta2) - 1, delta2;
                exp(i*beta3) - 1, delta3 ];
                
% Cramers rule for Ax=B using the 2 vector loop equations
W = det(numerator1) / det(denominator);
Z = det(numerator2) / det(denominator);

%Ground Pivot
A0 = PP1 - Z - W;
%moving pivot
A1 = PP1 - Z;

% print out all information in nice format
fprintf('PP1 = %f + %f *i \n', real(PP1), imag(PP1));
fprintf('PP2 = %f + %f *i \n', real(PP2), imag(PP2));
fprintf('PP3 = %f + %f *i \n', real(PP3), imag(PP3));

fprintf('alpha2 = %f (rad)\n', alpha2);
fprintf('alpha3 = %f (rad)\n', alpha3);

fprintf('W = %f + %f *i \n', real(W), imag(W));
fprintf('Z = %f + %f *i \n', real(Z), imag(Z));

fprintf('A0 = %f + %f *i \n', real(A0), imag(A0));
fprintf('A1 = %f + %f *i \n', real(A1), imag(A1));

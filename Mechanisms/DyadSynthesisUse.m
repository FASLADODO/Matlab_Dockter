%Problem 3 part C
%Rod Dockter
%October 2014

%specify precision positions and angle
PP1 = -279 + 152*i;
PP2 = -176 + 117*i;
PP3 = -105 + 36*i;
theta1 = 22.87;
theta2 = 352.75;
theta3 = 287.61;

%specify beta angles for Dyad A
beta2 = 300;
beta3 = 68;

%Compute dyad info using custom function based on code for prob 2
[A0, A1, WA, ZA] = DyadSynthesisFunc(PP1,PP2,PP3,theta1,theta2,theta3,beta2,beta3);

%specify beta angles for Dyad B
beta2 = 25;
beta3 = 57;

%Compute dyad info using custom function based on code for prob 2
[B0, B1, WB, ZB] = DyadSynthesisFunc(PP1,PP2,PP3,theta1,theta2,theta3,beta2,beta3);

%print out vals
fprintf('A0 = %f + %f *i \n', real(A0), imag(A0));
fprintf('A1 = %f + %f *i \n', real(A1), imag(A1));
fprintf('B0 = %f + %f *i \n', real(B0), imag(B0));
fprintf('B1 = %f + %f *i \n', real(B1), imag(B1));
fprintf('input W = %f + %f *i \n', real(WA), imag(WA));
fprintf('input Z = %f + %f *i \n', real(ZA), imag(ZA));
fprintf('output W = %f + %f *i \n', real(WB), imag(WB));
fprintf('output Z = %f + %f *i \n', real(ZB), imag(ZB));

% This yields these values
% A0 = -251.937948 + -42.893825 *i 
% A1 = -223.558398 + -10.210592 *i 
% B0 = -298.165343 + 185.315481 *i 
% B1 = -324.520955 + 51.833125 *i 
% input W = 28.379550 + 32.683232 *i 
% input Z = -55.441602 + 162.210592 *i 
% output W = -26.355611 + -133.482357 *i 
% output Z = 45.520955 + 100.166875 *i 
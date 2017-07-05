%%%%% Homework # 1, ME 8287

%% 2.8.3 a)

syms w1 w2 w3 theta

wHat = [0,-w3,w2;
    w3,0,-w1;
    -w2,w1,0];

disp('2.8.3 a) ')
[V,D] = eig(wHat)

%% 2.8.3 b)

R = expm(wHat*theta);

simpeig = subs(eig(R),(- w1^2 - w2^2 - w3^2)^(1/2),i); %exploiting ||W|| = 1 property

disp('2.8.3 b) 1')
simpeig = simplify(simpeig)

% getting eigen vectors for lambda = 1;

%[V,D] = eig(R) %returns no value

%instead by inspection we guess that the eigenvector for lambda = 1 is just
%the W unit vector 

W = [w1;w2;w3];

is_zero = (R - eye(3))* W;

is_zero = simplify(is_zero);

simpzero= subs(is_zero,(- w1^2 - w2^2 - w3^2)^(1/2),i);

simpzero = simplify(simpzero);

disp('2.8.3 b) 2')
pretty(simpzero) %Simplifies to zero hurray!

%% 2.8.10 b)

syms w theta real

wHat = [0, -w;
        w, 0];
    
disp('2.8.10 b)')
R = expm(wHat*theta)

% pretty(R) %simplifies the eulers

%% 2.8.10 c)

syms w theta real

R = [cos(theta*w),-sin(theta*w);
    sin(theta*w),cos(theta*w)];

wHat = [0, -w;
        w, 0];
    
result = R*wHat*transpose(R);

disp('2.8.10 c)')
simplify(result)
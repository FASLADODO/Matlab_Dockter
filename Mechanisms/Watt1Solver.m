% Watt 1 motion generation + prescribed timing solver
% Rod Dockter


% setup
clear all
deg2rad = (pi/180);
rad2deg = (180/pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify psi (angle between input positions) (deg)
psi2 = -90;
psi3 = -180;

% Precision positions real+i*im (x,y) (mm)
PP1 = 0 + 90*i;
PP2 = 25 + 92*i;
PP3 = 50 + 90*i;

% PP1 = 10 + 90*i;
% PP2 = 20 + 95*i;
% PP3 = 30 + 90*i;

% coupler angles
alpha2 = -0;
alpha3 = -0;
% alpha2 = -20;
% alpha3 = 20;

% ground pivot locations (measured from PP1) (mm)
PR_1 = 0 + 0*i;
PR_2 = 20 + 20*i;

% choose coupler lengths
Z3 = 25 + 20*i;
Z4 = -28 + 40*i;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%%%%

%convert to radians because matlab
psi2 = psi2 * deg2rad;
psi3 = psi3 * deg2rad;
alpha2 = alpha2 * deg2rad;
alpha3 = alpha3 * deg2rad;

% position differences (r*e^(i*theta) form )
delta = [ PP2 - PP1, PP3 - PP1 ];
delta_i = [ delta(1) - Z3*(exp(i*alpha2) - 1), delta(2) - Z3*(exp(i*alpha3) - 1) ];
delta_ii = [ delta(1) - Z4*(exp(i*alpha2) - 1), delta(2) - Z4*(exp(i*alpha3) - 1) ];

%ground pivot vectors dyad A
R1 = PP1 - Z3 - PR_1;
R2 = R1 + delta_i(1);
R3 = R1 + delta_i(2);


% Loop 1 Z1-Z2-Z3
% solve for phi (angle of Z2)
[phi2, phi3] = GroundSpecification(delta_i(1), delta_i(2), R1, R2, R3, psi2, psi3, 'Timing');
% Solve for Z1 and Z2
[Z1, Z2] = StandardDyadCramers(phi2, phi3, psi2, psi3, delta_i(1), delta_i(2));

%ground pivot vectors dyad B
R1 = PP1 - Z4 - PR_1;
R2 = R1 + delta_ii(1);
R3 = R1 + delta_ii(2);


% Loop 2 Z8-Z9-Z4
% solve for gamma (angle of Z9)
[gamma2, gamma3] = GroundSpecification(delta_ii(1), delta_ii(2), R1, R2, R3, psi2, psi3, 'Timing');
% Solve for Z8 and Z9
[Z8, Z9] = StandardDyadCramers(gamma2, gamma3, psi2, psi3, delta_ii(1), delta_ii(2));

%ground pivot vectors dyad B
R1 = PP1 - Z4 - PR_2;
R2 = R1 + delta_ii(1);
R3 = R1 + delta_ii(2);


% Loop 3 Z6-Z5-Z4
% solve for beta (angle of Z6)
[beta2, beta3] = GroundSpecification(delta_ii(1), delta_ii(2), R1, R2, R3, gamma2, gamma3, 'Motion');
% Solve for Z5 and Z6
[Z6, Z5] = StandardDyadCramers(gamma2, gamma3, beta2, beta3, delta_ii(1), delta_ii(2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Superfluous vectors
Z10 = Z1 - Z8;
Z11 = Z3 - Z4;
Z7 = Z5 - Z9;


% Plotting points
P0 = PR_1;
P1 = PR_2;
P2 = PR_1 + Z1;
P3 = PR_1 + Z8;
P4 = PR_2 + Z6;
P5 = PR_1 + Z1 + Z2;
P6 = PR_2 + Z6 + Z5;
P7 = PR_1 + Z1 + Z2 + Z3;


% Plot
figure(1)
plot(real(PP1),imag(PP1),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
hold on
plot(real(PP2),imag(PP2),'s','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
hold on
plot(real(PP3),imag(PP3),'d','MarkerSize',7,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(P0),imag(P0),'p','MarkerSize',12,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(P1),imag(P1),'p','MarkerSize',12,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(P2),imag(P2),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(P4),imag(P4),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(P3),imag(P3),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(P5),imag(P5),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(P6),imag(P6),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
fill( [real(P0), real(P2), real(P3) ], [imag(P0), imag(P2), imag(P3)] ,'g', 'FaceAlpha',0.6);
hold on
fill( [real(P3), real(P4), real(P6) ], [imag(P3), imag(P4), imag(P6)] ,'b', 'FaceAlpha',0.6);
hold on
fill( [real(P5), real(P6), real(P7) ], [imag(P5), imag(P6), imag(P7)] ,'c', 'FaceAlpha',0.6);
hold on
plot([real(P1), real(P4)],[imag(P1), imag(P4)],'-','LineWidth',2,'Color',[0 0 1]);
hold on
plot([real(P2), real(P5)],[imag(P2), imag(P5)],'-','LineWidth',2,'Color',[0 1 0]);

hold off

axis square equal
title('Watt 1')
xlabel('X (real)')
ylabel('Y (imaginary)')
legend('PP1','PP2','PP3','A0','B0','A1','B1','Location','BestOutside')


%%%%%%%%%%%%%%%%%%%%%%%%%%% Print stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%print out vals
fprintf('Mechanism Dimensions \n');
fprintf('A0 = %f + %f *i \n', real(P0), imag(P0));
fprintf('B0 = %f + %f *i \n', real(P1), imag(P1));
fprintf('Z1 = %f + %f *i \n', real(Z1), imag(Z1));
fprintf('Z2 = %f + %f *i \n', real(Z2), imag(Z2));
fprintf('Z3 = %f + %f *i \n', real(Z3), imag(Z3));
fprintf('Z4 = %f + %f *i \n', real(Z4), imag(Z4));
fprintf('Z5 = %f + %f *i \n', real(Z5), imag(Z5));
fprintf('Z6 = %f + %f *i \n', real(Z6), imag(Z6));
fprintf('Z7 = %f + %f *i \n', real(Z7), imag(Z7));
fprintf('Z8 = %f + %f *i \n', real(Z8), imag(Z8));
fprintf('Z9 = %f + %f *i \n', real(Z9), imag(Z9));
fprintf('Z10 = %f + %f *i \n', real(Z10), imag(Z10));
fprintf('Z11 = %f + %f *i \n', real(Z11), imag(Z11));
fprintf('alpha2 (rad)= %f \n', alpha2);
fprintf('alpha3 (rad)= %f \n', alpha3);
fprintf('beta2 (rad)= %f \n', beta2);
fprintf('beta3 (rad)= %f \n', beta3);
fprintf('psi2 (rad)= %f \n', psi2);
fprintf('psi3 (rad)= %f \n', psi3);
fprintf('phi2 (rad)= %f \n', phi2);
fprintf('phi3 (rad)= %f \n', phi3);
fprintf('gamma2 (rad)= %f \n', gamma2);
fprintf('gamma3 (rad)= %f \n', gamma3);

%% Print out positions


fprintf('P0 = %f + %f *i \n', real(P0), imag(P0));
fprintf('P1 = %f + %f *i \n', real(P1), imag(P1));
fprintf('P2 = %f + %f *i \n', real(P2), imag(P2));
fprintf('P3 = %f + %f *i \n', real(P3), imag(P3));
fprintf('P4 = %f + %f *i \n', real(P4), imag(P4));
fprintf('P5 = %f + %f *i \n', real(P5), imag(P5));
fprintf('P6 = %f + %f *i \n', real(P6), imag(P6));
fprintf('P7 = %f + %f *i \n', real(P7), imag(P7));






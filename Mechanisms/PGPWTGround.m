% Homework 4
% PGWPT with ground pivot specification
% Rod Dockter, Nov. 2014

% setup
clear all
deg2rad = (pi/180);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify Beta (angle between input positions) (deg)
beta2 = 36;
beta3 = 120;

% Precision positions real+i*im (x,y) (mm)
PP1 = -279 + 213*i;
PP2 = -88 + 163*i;
PP3 = -63 - 95*i;

% ground pivots (measured from PP1) (mm)
PR_1 = -495 - 167*i;
PR_2 = -257 - 257*i;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Implementation %%%%%%%%%%%%%%%%%%%%%%%%%
%Need to find ZA,WA,ZB,WB, alpha_j, and psi_j

% convert to radians because matlab
beta2 = beta2 *deg2rad;
beta3 = beta3 *deg2rad;

% position differences (r*e^(i*theta) form )
delta2 = PP2 - PP1;
delta3 = PP3 - PP1;


%%%%%%% start with dyad A, solve as ground pivot specification %%%%%%%%
%Dyad A is what I will call the dyad with prescribed input angle beta

% Start by finding alpha values using ground pivot methods
% This is a timing dyad ground specification 

%ground pivot vectors
R1 = PP1 - PR_1;
R2 = PP2 - PR_1;
R3 = PP3 - PR_1;


% Solving for alpha values using custom function
[alpha2, alpha3] = GroundSpecification(delta2, delta3, R1, R2, R3, beta2, beta3, 'Timing');


% now solve for dyad A (ZA and WA) with Cramers Rule Matrices
[WA, ZA] = StandardDyadCramers(alpha2, alpha3, beta2, beta3, delta2, delta3);




%%%%%% Now solve dyad B, solve as ground specification %%%%%
% input angle is not specified but now alpha angle is known

% Finding ground angle values for second dyad using ground pivot spec
% This is a motion dyad ground specification 


%ground pivot vectors
R1 = PP1 - PR_2;
R2 = PP2 - PR_2;
R3 = PP3 - PR_2;

% Solving for psi values using custom function
[psi2, psi3] = GroundSpecification(delta2, delta3, R1, R2, R3, alpha2, alpha3, 'Motion');


% now solve for dyad B (ZB and WB) with Cramers Rule Matrices
[WB, ZB] = StandardDyadCramers(alpha2, alpha3, psi2, psi3, delta2, delta3);


%moving pivots
A1 = PR_1 + WA;
B1 = PR_2 + WB;

%ground pivots (given)
A0 = PR_1;
B0 = PR_2;


%%%%%%%%%Check values %%%%%%%%%%%%%%%%%%%%%%%
% Moving pivots for 2 and 3 position
A2 = PR_1 + WA *exp(i*beta2);
A3 = PR_1 + WA *exp(i*beta3);

B2 = PR_2 + WB *exp(i*psi2);
B3 = PR_2 + WB *exp(i*psi3);

fprintf('check that W and Z values are same for all PP \n');

fprintf('Z1_1 = %f \n', abs(PP1 - A1));
fprintf('Z1_2 = %f \n', abs(PP2 - A2));
fprintf('Z1_3 = %f \n', abs(PP3 - A3));

fprintf('W1_1 = %f \n', abs(A1 - PR_1));
fprintf('W1_2 = %f \n', abs(A2 - PR_1));
fprintf('W1_3 = %f \n', abs(A3 - PR_1));

fprintf('Z2_1 = %f \n', abs(PP1 - B1));
fprintf('Z2_2 = %f \n', abs(PP2 - B2));
fprintf('Z2_3 = %f \n', abs(PP3 - B3));

fprintf('W2_1 = %f \n', abs(B1 - PR_2));
fprintf('W2_2 = %f \n', abs(B2 - PR_2));
fprintf('W2_3 = %f \n', abs(B3 - PR_2));
fprintf('\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%print out vals
fprintf('Mechanism Dimensions \n');
fprintf('A0 = %f + %f *i \n', real(A0), imag(A0));
fprintf('A1 = %f + %f *i \n', real(A1), imag(A1));
fprintf('B0 = %f + %f *i \n', real(B0), imag(B0));
fprintf('B1 = %f + %f *i \n', real(B1), imag(B1));
fprintf('DyadA W = %f \n', abs(WA));
fprintf('DyadA Z = %f \n', abs(ZA));
fprintf('DyadB W = %f \n', abs(WB));
fprintf('DyadB Z = %f \n', abs(ZB));
fprintf('alpha2 (rad)= %f \n', alpha2);
fprintf('alpha3 (rad)= %f \n', alpha3);
fprintf('beta2 (rad)= %f \n', beta2);
fprintf('beta3 (rad)= %f \n', beta3);
fprintf('psi2 (rad)= %f \n', psi2);
fprintf('psi3 (rad)= %f \n', psi3);


%%%%% Plotting %%%%%%
figure(1)
plot(real(PP1),imag(PP1),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
hold on
plot(real(PP2),imag(PP2),'s','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
hold on
plot(real(PP3),imag(PP3),'d','MarkerSize',7,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(A0),imag(A0),'p','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(B0),imag(B0),'p','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(A1),imag(A1),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(B1),imag(B1),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot([real(A0), real(A1)],[imag(A0), imag(A1)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(B0), real(B1)],[imag(B0), imag(B1)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(B1), real(PP1)],[imag(B1), imag(PP1)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(A1), real(PP1)],[imag(A1), imag(PP1)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(A1), real(B1)],[imag(A1), imag(B1)],'-','LineWidth',2,'Color',[0 0 0]);

hold off

axis equal
title('PGWPT + Ground Spec (PP1)')
xlabel('X (real)')
ylabel('Y (imaginary)')
legend('PP1','PP2','PP3','A0','B0','A1','B1','Location','BestOutside')

figure(2)
plot(real(PP1),imag(PP1),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
hold on
plot(real(PP2),imag(PP2),'s','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
hold on
plot(real(PP3),imag(PP3),'d','MarkerSize',7,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(A0),imag(A0),'p','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(B0),imag(B0),'p','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(A2),imag(A2),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(B2),imag(B2),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot([real(A0), real(A2)],[imag(A0), imag(A2)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(B0), real(B2)],[imag(B0), imag(B2)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(B2), real(PP2)],[imag(B2), imag(PP2)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(A2), real(PP2)],[imag(A2), imag(PP2)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(A2), real(B2)],[imag(A2), imag(B2)],'-','LineWidth',2,'Color',[0 0 0]);
hold off

axis equal
title('PGWPT + Ground Spec (PP2)')
xlabel('X (real)')
ylabel('Y (imaginary)')
legend('PP1','PP2','PP3','A0','B0','A1','B1','Location','BestOutside')

figure(3)
plot(real(PP1),imag(PP1),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
hold on
plot(real(PP2),imag(PP2),'s','MarkerSize',7,'MarkerFaceColor',[1 0 1],'Color',[1 0 1]);
hold on
plot(real(PP3),imag(PP3),'d','MarkerSize',7,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(A0),imag(A0),'p','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(B0),imag(B0),'p','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot(real(A3),imag(A3),'o','MarkerSize',8,'MarkerFaceColor',[0 1 0],'Color',[0 1 0]);
hold on
plot(real(B3),imag(B3),'o','MarkerSize',8,'MarkerFaceColor',[0 0 1],'Color',[0 0 1]);
hold on
plot([real(A0), real(A3)],[imag(A0), imag(A3)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(B0), real(B3)],[imag(B0), imag(B3)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(B3), real(PP3)],[imag(B3), imag(PP3)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(A3), real(PP3)],[imag(A3), imag(PP3)],'-','LineWidth',2,'Color',[0 0 0]);
hold on
plot([real(A3), real(B3)],[imag(A3), imag(B3)],'-','LineWidth',2,'Color',[0 0 0]);
hold off

axis equal
title('PGWPT + Ground Spec (PP3)')
xlabel('X (real)')
ylabel('Y (imaginary)')
legend('PP1','PP2','PP3','A0','B0','A1','B1','Location','BestOutside')






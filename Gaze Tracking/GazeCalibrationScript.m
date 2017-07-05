%Script for doing gaze calibration
%Import workspace first

gazedata1 = dlmread('gazedata_2.txt'); %change to proper file, delete header row

kin_data = gazedata1(:,1:13);

[tipL,tipR] = EDGE_Kinematics(kin_data,2,0);

XR = tipR(:,1);
YR = tipR(:,2);
ZR = tipR(:,3);

XL = tipL(:,1);
YL = tipL(:,2);
ZL = tipL(:,3);

Gx = gazedata1(:,14);
Gy = gazedata1(:,15);

scatter3(XR,YR,ZR);

%scatter(Gx,Gy);

%%

% using the following form phi = (D^T*D)^-1*D^T * G
%D is 4xn data matrix

%Columns of transformation matrix
phi_1 = [0;0;0;0];
phi_2 = [0;0;0;0];
phi_3 = [0;0;0;0];
phi_4 = [0;0;0;0];

%Data matrix
D = [XR,YR,ZR,ones(length(XR),1)];

G_1 = Gx;
G_2 = Gy;
G_3 = zeros(length(Gx),1);
G_4 = ones(length(Gx),1);

%Compute the transformation matrix

premultiply = inv(transpose(D)*D)*transpose(D);

phi_1 = premultiply*G_1;
phi_2 = premultiply*G_2;
phi_3 = premultiply*G_3;
phi_4 = premultiply*G_4;

trans_matrix = [transpose(phi_1);transpose(phi_2);transpose(phi_3);transpose(phi_4)]


%%

%calibrated data

%mat(row,column)

data_out = zeros(4,length(XR));

for i = 1:length(XR)
    data_out(:,i) = trans_matrix*transpose(D(i,:));
end

%%

%plotting

scatter(Gx,Gy)


hold on

scatter(data_out(1,:),data_out(2,:))

hold off
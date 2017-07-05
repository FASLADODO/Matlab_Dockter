%function icpTest
% t=[-1:.001:4]';
% tfew = [0:.5:3];
% D = pCurve(t);
% Dfew = pCurve(tfew);
clear global;
corners =[ 0 0 0 ; 1 0 0; 1 5 0 ; 0 7 0; 0 7 3];
m=mean(corners);
r= max(range(corners))*[-.5 .5];
ax3d = [m(1)+r(1) m(1)+r(2) m(2)+r(1) m(2)+r(2) m(3)+r(1) m(3)+r(2)];

D = refBlockPointGen(corners, 10);
Dfew = corners; %refBlockPointGen(corners, 1);
np = 40; % number of points to skip
mp = round(size(D,1)/2)-round(np/2); 
m1 = 1:mp; 
m2 = (mp+np):size(D,1);

% crap =union(find(isnan(D(:,1))), find(isnan(D(:,2))));
% D(crap,:) =[];
% t(crap) = [];

% Add noise to test data
Dn = D +randn(size(D))*.1;
Dnpartial = Dn([m1 m2:end ],:);
Dn=Dnpartial;
%D = Dn([m1 m2:end ],:);

% Rotate and offset test data to generate "samples"
r = pi/3;
p0 = [5 5 0];
R =angle2dcm(r, 0, 0);
T = [R p0'];
Ds = [T*[Dn ones( length(Dn),1) ]']' ;

%% Run ICP algorithm
%[Rfit, Tfit] = icp(model,data)
%   OUTPUT:
%  
%     R - rotation matrix and
%     T - translation vector accordingly so
%  
%             newdata = R*data + T .
%[Rfit, Tfit] = icp(D,Ds);
max_iter=400;
min_iter=40;
fitting=2; %[2,w]
thres=1e-5;
init_flag=1;
tes_flag=1;

tic;
[Rfit, Tfit] = icp(D, Ds,...
    max_iter,min_iter,fitting,thres,init_flag,tes_flag);
T2 = Rfit'*Tfit ;% In original coordinates
D2=zeros(size(Ds));
for i=1:length(Ds)
    D2(i,:) = [(Rfit)*Ds(i,:)' + Tfit]';
end
toc
%Dfit = inv(Rfit)*[(Ds) zeros(size(Ds,1),1)]';
%Dfit=Dfit';





%% Now Plot stuff ...
figure(10); clf
%plot3(t, D(:,1), D(:,2)); hold on
%plot3(D(:,1), D(:,2),D(:,3),'o'); hold on
grid on; 
axis square ;axis(ax3d)

%figure(20);clf
plot3(Dfew(:,1), Dfew(:,2),Dfew(:,3),'k-','linewidth',1 ); % Shape
grid on; hold on; axis equal;
plot3(Dfew(:,1), Dfew(:,2),Dfew(:,3),'k.','markersize',20 ); % Control Points
plot3(D(:,1), D(:,2),D(:,3),'k>','markersize', 4); % Ref Points
axis square ;axis(ax3d)



%plot(Dn(:,1), Dn(:,2),'.')
plot3(Ds(:,1), Ds(:,2),Ds(:,3),'b.'); % Sampled data
drawnow;



plot3(D2(:,1), D2(:,2),D2(:,3),'ro', 'markersize', 3); % fitted data
title('\bfIterative Closest Point (ICP) Method');
xlabel('x [cm]')
ylabel('y [cm]')
legend('Shape','Control Point', 'Ref point', 'Samples', 'Fit', 'location', 'best')
axis square ;axis(ax3d)

disp('Original offset, guessed offset, diff:')
disp([p0' -T2 T2+p0'])
disp('Original R, guessed R:')
disp([R' inv(Rfit) ])

return
%% %% define parametric curve:
% function D = pCurve(t)
% corners =[ 0 0; 1 0; 1 5; 0 7];
% %segments of the parameter t in p(t); 
% % first and last, between [-Inf tSeg(1)] and finally: [tSeg(end) Inf], get
% % mapped to [];  
% % then [tSeg(k) tSeg(k+1)] get mapped to corners' segments linearly
% tSeg = [0 1 2 3.00001 ]; 
% D=zeros(length(t),2)*NaN;
% %t(find(t<tSeg(1)))=[];
% %t(find(t>tSeg(end)))=[];
% 
% for i=[2:(size(corners,1))];
%     ts = intersect( find(t>=tSeg(i-1)), find(t<=tSeg(i)));
%     for tt=ts;
%         
%         D(tt, :)= [ ...t(tt) tt...
%             (corners(i,1)-corners(i-1,1))*(t(tt)-tSeg(i-1))./(tSeg(i)-tSeg(i-1))+corners(i-1,1)...
%             (corners(i,2)-corners(i-1,2))*(t(tt)-tSeg(i-1))./(tSeg(i)-tSeg(i-1))+corners(i-1,2)...
%             ];
%     end
%end
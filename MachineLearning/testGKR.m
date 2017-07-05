% Taken from here:
% http://youngmok.com/gaussian-kernel-regression-with-matlab-code/

% test gaussian kernel regression
% functionized version: gaussian_kern_reg.m

%% 1-D case


x = 1:100; % training data x
y = sin(x/10)+(x/50).^2+0.2*randn(1,100);   % training data y
h=5; % kernel bandwidth

for i=1:100
    xs(i)=i;
    %     ys(i)=gaussian_kern_reg(xs(i),x,y,h);
    
    % Gaussian kernel function
    z = (xs(i)- x ) /h ;
    k = exp(-z.*z/2)/sqrt(2*pi);
    
    ys(i,:) = sum(k.*y)/sum(k);
end

figure;hold on; 
plot(x,y,'.');
plot(xs,ys,'r-');

%% 2-D case

nn = 100;

%Dataset Generation
x = rand(nn,2)-0.5;
y = (x(:,1)-0.5).^2 + x(:,2)+  0.1*rand(nn,1);

% Gaussian Kernel Bandwidth Setting
h=[0.1;0.1]; 

xl1 = linspace(min(x(:,1)),max(x(:,1)),50);
xl2 = linspace(min(x(:,2)),max(x(:,2)),50);

%Prediction point
[xx1, xx2] = meshgrid(xl1,xl2);

for i=1:size(xx1,1)
    for j=1:size(xx1,2)
        xs=[xx1(i,j),xx2(i,j)];
        ys(i,j)=gaussian_kern_reg(xs,x,y,h); % Prediction
    end
end

% Result Plot
figure;hold on; 
plot3(x(:,1),x(:,2),y,'.')  % training dataset plot
mesh(xx1,xx2,ys)            % prediction point plot
legend('Training data','Predicted Surface')
grid ; view(3);


%% 2-D weird case


nn = 100;
scale = 1.5;

x1 = [1:nn]' * scale; % training data x
x2 = [1:nn]'; % training data x

%make a grid
xlin = linspace(min(x1),max(x1),50);
ylin = linspace(min(x2),max(x2),50);
[X,Y] = meshgrid(xlin,ylin);

x = [X(:), Y(:)];

zy = sin(x(:,1)./10) + sin(x(:,2)./2) + (x(:,1)./50).^2 + 0.2*randn(length(x),1);   % training data y
f = scatteredInterpolant(x(:,1),x(:,2),zy,'linear', 'none');
Z = f(X,Y);

h=2; % kernel bandwidth

mesh(X,Y,Z);

%%
[NN,SS] = size(x);

for i=1:NN
    xs(i,:)=x(i,:);
    ys(i,:)=gaussian_kern_reg(xs(i,:),x,zy,h); % Prediction
    
    % Gaussian kernel function
%     z = NormRowWise( ( repmat(xs(i,:),NN,1) - x ) /h );
%     k = exp(-z.*z/2)/sqrt(2*pi);
%     
%     ys(i,:) = sum(k.*zy)/sum(k);

end

figure;hold on; 
h=mesh(X,Y,Z);
set(h,'facecolor','none')
scatter3(xs(:,1),xs(:,2),ys,'r.');
hold off

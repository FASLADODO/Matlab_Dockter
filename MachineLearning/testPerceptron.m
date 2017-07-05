%%perceptron

load fisheriris

DataM = meas;

[YM, classes] = NumericClassLabels(species);
newclass = [-1,1,10];

for ii = 1:length(classes)
    
    datac{ii} = DataM(YM == ii,:);
    temp = YM(YM  == ii);
    yc{ii} = ones(length(temp),1)*newclass(ii);
end


figure
h1 = scatter3(datac{1}(:,1),datac{1}(:,2),datac{1}(:,3),'r.');
hold on
h2 = scatter3(datac{2}(:,1),datac{2}(:,2),datac{2}(:,3),'b.');
hold on
h3 = scatter3(datac{3}(:,1),datac{3}(:,2),datac{3}(:,3),'g.');
hold off
legend([h1(1),h2(1),h3(1)],classes)

%%

Data = [];
Data = [datac{1};datac{2}]; %just tw oclasses for now
Data = [ones(length(Data),1),Data]; %add in a ones column as bias (b) term
dout = [yc{1} ; yc{2} ]; %perceptron value -1 or 1

maxiter = 1000;
gamma = 0.005; %small enoug rms
alpha = 0.0001; %learning rate

Weights = perceptron(Data, dout, alpha, maxiter,gamma);

Weights

dcheck = Data*Weights;

figure
plot(dcheck)

ClassCheck = zeros(length(dout),1);

ClassCheck(dcheck < 0) = 1;
ClassCheck(dcheck > 0) = 2;


%% Try linear

nn=100;
x0 = linspace(1,10,nn)';
y = 2.5*x0 + 3*ones(nn,1) + rand(nn,1)*0.1;

scatter(x0,y)

Data = [ones(nn,1),x0];

maxiter = 1000;
gamma = 0.0005; %small enoug rms
alpha = 0.001; %learning rate

Weights = perceptron(Data, y, alpha, maxiter,gamma)

lsw = pinv(Data)*y



nlinfunc = @(Param,X) (Param(2)*X(:,2) + exp(Param(1)*X(:,1)) ) ; 


params1 = [2,3];
params2 = [3,4];
paramact = [params1; params2];

samples = 10;
tt = 50;
noise = 0.01;
pnoise = 0.4;

Data = [];
for gg = 1:2
    for ss = 1:samples
       theta = [1:tt]' + randn(tt,1)*noise;  
       thetadot = diff(theta);
       dtemp = [theta, thetadot];
       
       ptemp = paramact(gg,:);
       
       
    end

end
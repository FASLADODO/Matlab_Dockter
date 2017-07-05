function Weights = perceptron(Data, Y, alpha, maxiter,gamma)
%Data is column vectors with a ones vector for bias
%Y = -1 or 1
%alpha is learning rate 0 < alpha < 1
%maxiter is maximum number of iterations
%gamma is threshold below which RMS will exit iterations (0.001)

dout = Y;

[NN,SS] = size(Data);

lambda = 0.001; %params not changing

W = rand(SS,1)*4; %rando
prevW = W;
for tt = 1:maxiter
    for ii = 1:NN
        xj = Data(ii,:); %current data sample
        f_z = xj * W;
        res = dout(ii,1) - f_z;
        W = W + alpha.*res.*xj';
        
        if(norm(prevW - W) < gamma)
           %fprintf('\r We not movin much \n')
        end
        
        prevW = W;
    end
    
    residual = (1/NN) *mean( abs( dout - (Data*W) ) );
    
    if(residual < gamma)
       break; 
    end
    
end

fprintf('Total Iterations: %d, Residual: %f \n',tt, residual)

Weights = W; %return em
function [mua,mub,sigabbbinv,siga_b] = ConditionalGaussian(xa,xb)
%See page Page 87 Bishop
%Returns the gaussian for parameters for P(xa|xb)
%xa and xb must both be column vectors with each column as a dimension

%compute mua_b for vector of xb as:
%mua_b = repmat(mua',1,NNb) + sigab*inv(sigbb)*(xb - repmat(mub,NNb,1))';
%mua_b = mua_b';

    %FUNCTION START
    [NNa,SSa] = size(xa);
    [NNb,SSb] = size(xb);

    %marginal distributions
    xall = [xa,xb];

    mua = mean(xa); %mean of group a
    mub = mean(xb); %mean of group b

    sigall = cov(xall)
    muall = mean(xall)

    %break down sigall into parts (pg 87 bishop)
    sigaa = sigall(1:SSa,1:SSa);
    sigbb = sigall(SSa+1:end,SSa+1:end);
    sigab = sigall(1:SSa,SSa+1:end);
    sigba = sigall(SSa+1:end,1:SSa); %=sigab'

    %get inverses of covariance parts
    AL11 = inv(sigaa - sigab*inv(sigbb)*sigba);
    AL12 = -inv(sigaa - sigab*inv(sigbb)*sigba)*sigab*inv(sigbb);


    %conditional parameters
    sigabbbinv = sigab*inv(sigbb); %for online use
    siga_b = sigaa - sigab*inv(sigbb)*sigba;
end
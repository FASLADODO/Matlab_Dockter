function [pa_b] = ConditionalGaussianDist(Xb,mua,mub,sigabbbinv,siga_b)
%See page Page 87 Bishop
%sigabbbinv = sigaa*inv(sigbb) to speed things up

    [NN,SS] = size(Xb);

    mua_b = repmat(mua',1,NN) + sigabbbinv*(Xb - repmat(mub,NN,1))';
    mua_b = mua_b';
    
    sigin = inv(siga_b);
    ss = length(siga_b);

    P = (1/ ( (2*pi)^(ss/2) .*sqrt(norm(siga_b)) ) ).* exp(-(1/2).* sum((Xb - mua_b)' .* (sigin*(Xb - mua_b)'), 1));
    
    pa_b = P';
end
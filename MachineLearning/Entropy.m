function H = Entropy(P)
%https://en.wikipedia.org/wiki/Entropy_(information_theory)
%If there are multiple columns, do entropy seperately for each column and
%sum

    [~,SS] = size(P);

    H = 0;
    for ii = 1:SS
        ptemp = P(:,ii);
        if(ptemp == 0)
            H = H + 0;
        else
            H = H + -sum(ptemp.*log2(ptemp));
        end
    end
end
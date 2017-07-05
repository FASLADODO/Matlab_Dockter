function H = EntropyVector(P)
%https://en.wikipedia.org/wiki/Entropy_(information_theory)
%fora binary variable

    H = -P.*log2(P)  - (1-P).*log2(1-P);

end
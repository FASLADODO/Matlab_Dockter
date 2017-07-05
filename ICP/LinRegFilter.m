N=4;

Arr = linspace(1,N,N);
vals = zeros(1,N);     

B = 1/(1-(1/N)*sum(Arr)+N);

M = (1-N*B)/sum(Arr);

for i = 1:N
    vals(i) = M*Arr(i) + B;
end

B
Arr
vals


% filter(b(y coefficients), a(x coefficients), DATA)
freqz(0, vals, 0:1000) 
%%

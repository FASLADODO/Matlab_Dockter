%root locus

s = tf('s');

y = ((s+3))/((s+2)*(s+3)*(s))

rlocus(y)

%%

w = -1 + 1*i;

G = (w+3)/(w*(w+3)*(w+2));

K = 1/abs(G)

%%

syms s

expand((s+1-i)*(s+1+i))

%% jw crossing

syms w k

solve(-k*w^2+4*k+6*w^4-60*w^2==0,4*k*w+w^5-35*w^3+250*w == 0)

%%

syms s

solve(s^2-2*s+2 == 0)
function[pi]=make_pi(N)
hits=0;
for j=1:N
    X=rand;
    Y=rand;
    R=X.^2+Y.^2;
    inside=R<=1.;
    hits=hits+sum(inside);
end
ratio=hits/N;
pi=ratio*4;



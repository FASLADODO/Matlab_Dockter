%fermats last theorem

n = 5;
M = 100;


%get as and bs
a = 1:M;
b = 1:M;

%put them in cross wise grid
[ag,bg] = meshgrid(a,b);


cg = (ag.^n + bg.^n).^(1/n);

isint = cg - floor(cg);

mesh(ag,bg,isint)

min(min(isint))
max(max(isint))
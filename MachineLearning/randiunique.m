function r = randiunique(imax,m)
%get unique random integers
    r=randperm(imax);
    r=r(1:m);
end
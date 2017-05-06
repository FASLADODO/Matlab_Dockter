A = [1,2,3,4,5,6,7;
    12,523,6363,685,9673,234,213;
    8,8,8,8,83,36,89];

[row,col] = size(A);

parfor i = 1:row
   mm = mean(A(i,:)) 
   fprintf('mean = %f, row = %d',mm,i);
    
end
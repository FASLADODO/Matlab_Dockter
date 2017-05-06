function A = bubblesort( A )
% http://www.algorithmist.com/index.php/Bubble_sort
%stupid slow

N = length(A);

for i = 1 : N
    for j = 1:N - 1
       if A(j) > A(j+1)
          %swap em
          A([j,j+1]) = A([j+1,j]);
       end
    end
end

    
end
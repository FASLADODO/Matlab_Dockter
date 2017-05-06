function A = Partition(A)

     
     asize = length(A);
     ri = randi([1 asize],1,1);
     pivot = A(ri);

     left = [];
     right = [];
     for i = 1:asize
         if (A(i) < pivot)
             left = [left, A(i)];
         else
             right = [right, A(i)];
         end
     end
     
     A = [left,pivot,right];
end
function A = Quicksort(A)
%recursive implementation
%http://www.mathcs.emory.edu/~cheung/Courses/171/Syllabus/7-Sort/quick-sort1.html
%https://en.wikipedia.org/wiki/Quicksort

%fastest if data is not well behaved
    
    %get size and check if we're done with that guy
    asize = length(A);
    if(asize <= 1)
        return;
    end
    
    %get a random pivot
    ri = randi([1 asize],1,1);
    pivot = A(ri);

    %partition (YAY matlab)
    left = A(A < pivot);
    right = A(A >= pivot);

    %now sort each partition seperately (recursion is scary)
    left = Quicksort(left);
    right = Quicksort(right);

    %concatenate
    A = [left,right];

end
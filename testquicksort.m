%%try quicksort sorting algorithms

nn = 100;

%Data 
Data = randn(nn,1);

%time it
tic
A = Quicksort(Data);
toc

%compare to bubble (slow)
tic
B = bubblesort( Data );
toc

%compare to insertion sort
tic
[X_Sort,IDX_Sort] = InsertionSort(Data);
toc

figure
bar(1:nn,Data);

figure
bar(1:nn,A);

figure
bar(1:nn,B);
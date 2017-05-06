function A = Swap(A,id1,id2)
%simple array swap
    A([id1,id2]) = A([id2,id1]);
end
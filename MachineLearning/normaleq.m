function theta = normaleq(X, Y)
    %page 11 cs229-notes1
    theta = inv(X'*X)*X'*Y;
end
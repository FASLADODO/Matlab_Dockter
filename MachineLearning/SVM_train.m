function [W,B,SV] = SVM_train(X,Y)
%Solve SVM in matlab: http://www.csd.uwo.ca/~olga/Courses/CS434a_541a/Lecture11.pdf
    [NN,SS] = size(X);
    H = (X*X').*(Y*Y');
    f = -ones(NN,1);
    A = -eye(NN);
    a = zeros(NN,1);
    B = [Y' ; zeros(NN-1,NN) ];
    b = zeros(NN,1);

    marginz = 0.001;
    alpha = quadprog(H + eye(NN)*marginz, f, A, a, B, b);
    %Quadratic programming expects:
    % L = 0.5*alpha'*H*alpha + f'*alpha
    % A*alpha < a and B*alpha = b

    %floor noise
    thresher = 0.001;
    alpha(find(alpha < thresher)) = 0;
    %get support vectors
    idSV = find(alpha > 0);
    SV = X(idSV,:); %Supports vectors

    %find W 
    W = (alpha.*Y)'*X;
    B = mean((1./Y(idSV)) - X(idSV,:)*W'); %average of all W_0

end
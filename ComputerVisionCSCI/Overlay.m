function Combine = Overlay(Matrix1,Matrix2,overlap)
%%Adds two matrices together with a specified overlap to create larger
%%matrix
%%give square matrices PLZ
[M1,N1]=size(Matrix1);
[M2,N2]=size(Matrix2);
if M1 ~= N1 || M2 ~= N2
    error('Overlay:DimChk', 'Both matrices must be square.');
end
if M1 == N1 || M2 == N2
    newsize=M1+M2-overlap;
    point=M1-overlap+1;
    %%dummy matrix
    C = zeros(newsize);
    C(1:M1,1:M1) = Matrix1;
    C(point:newsize,point:newsize) = C(point:newsize,point:newsize) + Matrix2;
end
Combine=C;
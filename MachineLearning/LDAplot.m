function h = LDAplot(Data,W,C)
%W and C come from LDASimple() (Rods Weird LDA)

%C = MdlLinear.Coeffs(1, 2).Const
%W = MdlLinear.Coeffs(1, 2).Linear

f = @(x,y) C(2) + [x, y]*W(2,:)';

h = ezplot(f, [min(Data(:,1)) max(Data(:,1)) min(Data(:,2)) max(Data(:,2))]);

end
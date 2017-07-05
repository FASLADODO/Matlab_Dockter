function [xmean,ymean,xstd,ystd,rstd] = MCSIM(N)
A=zeros();
B=zeros();
ND=25.0;
dstd=0.0001;
sstd1=0.002;
sstd2=0.008;
sstd3=0.009;
sstd4=0.009;
% d=(randn(1,N)*dstd)+79.2;
% Y1=(randn(1,N)*sstd1)+-147.90199457749;
% Y2=(randn(1,N)*sstd2)+-107.8997683037364;
% Y3=(randn(1,N)*sstd3)+107.8997683037364;
% Y4=(randn(1,N)*sstd4)+147.90199457749;
d=(randn(1,N)*dstd)+79.544;
Y1=(randn(1,N)*sstd1)+15.16332;
Y2=(randn(1,N)*sstd2)+58.23686;
Y3=(randn(1,N)*sstd3)+265.8588;
Y4=(randn(1,N)*sstd4)+309.6434;
for j=1:N
    Y=[Y1(j), Y2(j),Y3(j),Y4(j)];
    X=[ND,ND + d(j),ND + d(j),ND];
    xbar = (1/numel(X))*sum(X);
    ybar = (1/numel(Y))*sum(Y);
    u = X-xbar;
    v = Y-ybar;
    u2 = sum(u.*u);
    uv = sum(u.*v);
    u3 = sum(u.*u.*u);
    uv2 = sum(u.*v.*v);
    v2 = sum(v.*v);
    v3 = sum(v.*v.*v);
    vu2 = sum(v.*u.*u);
    syms uc vc;
    f(1) = uc*u2+vc*uv-(1/2)*(u3+uv2);
    f(2) = vc*v2+uc*uv-(1/2)*(v3+vu2);
    [uc,vc] = solve(f(1), f(2));
    xc=vpa(uc+xbar);
    yc=vpa(vc+ybar);
    A(j)=xc;
    B(j)=yc;
end
A;
B;
figure(1)
title('X-center')
xlabel('X coordinate'), ylabel('frequency')
hist(A,50)
figure(2)
title('Y-center')
xlabel('Y coordinate'), ylabel('frequency')
hist(B,50)
xmean=mean(A);
ymean=mean(B);
xstd=std(A);
ystd=std(B);
rstd=sqrt(xstd^2+ystd^2);






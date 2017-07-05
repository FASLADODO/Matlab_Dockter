function [V,ss]=nystrom_ncut_2b(A,B);
% [V,ss]=nystrom_ncut_2b(A,B);

Nsamples = size(A,1);
Nothers = size(B,2);
Npix = Nsamples + Nothers;

% compute total connection weight 
disp('computing total connection weight...');
d1 = sum([A;B'],1);
d2 = sum(B,1) + sum(B',1)*pinv(A)*B;
d = [d1 d2]';
disp('done.')

% normalize
disp('normalizing...');
v = sqrt(1./d);
A = A.*(v(1:Nsamples)*v(1:Nsamples)');
B = B.*(v(1:Nsamples)*v(Nsamples+(1:Nothers))');
disp('done.')

% find eigenvectors via PCA/nystrom trick
disp('computing eigenvectors...')
[U,S,junk]=svd(A);
Asi=U*pinv(sqrt(S))*U';
Q=A+Asi*B*B'*Asi;
[U,L,junk]=svd(Q);
Va=[A;B']*Asi*U*pinv(sqrt(L));

for i = 1:Nsamples-1
  V(:,i) = Va(:,i+1)./Va(:,1);
end;
disp('done.')

ss=1-diag(L); % convert to (D-W)x=lambda Dx form
ss(1)=[]; % drop first eigenvalue (which is 1)


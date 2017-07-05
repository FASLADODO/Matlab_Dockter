%make random data
nn = 50;

sig1 = 2.0;
mu1 = 1;
d1 = randn(nn, 1)*sig1 + mu1;

sig2 = 0.2;
mu2 = 8;
d2 = randn(nn, 1)*sig2 + mu2;

%give it labels
G = [d1;d2];
L = [ones(nn,1)*1;ones(nn,1)*2];
DIRL = [1;2];

axvs = ones(2*nn,1);

%% 

%compute std
s1 = std(d1)
s2 = std(d2)


% if(s1 > s2)
%     grad = (s1./s2) .* sign(s2-s1)
% else
% 	grad = (s2./s1) .* sign(s1-s2)
% end
% threshstep = (atan(0.5*grad) / pi) + 0.5


grad = ((min([s1,s2])/max([s1,s2])) - 1)
threshstep = 0.5 + sign(s2-s1)*grad*0.08

% all of these suck
% threshstep = (3 ./ ((s2-s1))) + 0.5
% threshstep = s1 ./ (s1 + s2)
% threshstep = s2 ./ (2*s1)
% threshstep = s1 ./ (s2 + 1)
% threshstep = ((s1-s2) ./ (s1 + s2))
% threshstep = abs(s1-s2) ./ (s1 + s2)
% threshstep = 0.5

% sort and try thresholds
GS = sort(G);
GD = diff(GS);

accall = 0;
thrshB = 0;
bordz = [];
for ii = 1:length(GD)
   thrsh = GS(ii) + GD(ii)*threshstep;
   abv = G > thrsh;
   classG = DIRL(abv+1);
   corr = classG == L;
   accG = mean(corr);
   if(accG > accall)  %save best thresh
      accall = accG;
      thrshB = thrsh;
      bordz = [GS(ii), GS(ii+1)];
   end
end
thrshB

%plot it
figure
gscatter(axvs,G,L)
hold on
scatter(1,thrshB,'gx')
hold on
scatter(1,bordz(1),'ko')
hold on
scatter(1,bordz(2),'ko')
hold off
str = sprintf('threshstep = %f',threshstep);
title(str)

%% now it's weaponized

W = ones(2*nn,1)/(2*nn);
[Thresh,bestscore,IG] = DecisionStumpBasic(G,W,L)




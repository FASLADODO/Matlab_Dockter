% arcLength(d, th) computes the total pythagorean distance across the rows 
% of data d and uses th as a threshold, point-to-point distances less than 
% threshold do not contribute to computed arc length;
%  
%   Input:  d: (mxn) data where m= number of samples, n= dimension
%          th: (1x1) length delta threshold (min path segment legnth) OR
%             |(1xn) dim-specific delta threshold  
%
%   Output: L (mx1) array of lengths, starting at 0 for d(1,:);
function L = arcLength(d,th)

df = diff(d);

% zero out diffs that fall below dim-specific thresholds...th is 1xn
if length(th)~=1
    df(abs(df)<=repmat(th, size(df,1), 1)) = 0;   
end

% compute arc Length
aL = sqrt(sum(df.^2 ,2));
%aLo=aL;

% zero out any path lengths below th of size 1x1
if length(th)==1
    aL(abs(aL)<=th) = 0;
end

% disp([aLo aL]);

% compute cumulative sum and pad with zero in beginning. 
L = cumsum( [0 ; aL],1);
function [centers,assignments,error] = kmwrapper(data,k,n,t,count);
er=inf;

for c=1:count
  c
  [cen,a,e]=km(data,k,n,t);
  if sum(e) < er;
    er=sum(e)
    centers=cen;
    assignments=a;
    error=e;
  end
end

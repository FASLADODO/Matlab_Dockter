function PlotActuators(starter,vector)

[dim,kk] = size(starter);
[dimc,kkc] = size(vector);

if(dim ~= dimc || kk ~= kkc)
    error('vectors must be same size')
end

ender = starter + vector;

for ii = 1:kk
    plot3([starter(1,ii),ender(1,ii)],[starter(2,ii),ender(2,ii)],[starter(3,ii),ender(3,ii)],'k-','LineWidth',5)
    hold on
end

end
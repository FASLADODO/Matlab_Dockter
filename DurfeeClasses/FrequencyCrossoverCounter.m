
crosses=0;
maxamp = 0;
minamp = 0;

avg = mean(vals);
for j = 1:length(vals) - 1
    old = vals(j) - avg;
    new = vals(j+1) - avg;
    if sign(old) ~= sign(new)
        crosses = crosses + 1;
    end
    if old > maxamp
        maxamp = old;
    end
    if old < minamp
        minamp = old;
    end
end

avg
crosses
minamp
maxamp


%%
crosses=0;
maxamp = 0;
minamp = 0;
flipflop = 0;

threshold = 100;

avg = mean(vals);
for j = 1:length(vals)
    data = vals(j) - avg;
    
    if flipflop == 0
        if data > threshold
            crosses = crosses + 1;
            flipflop = -1;
        end
        if data <  -threshold
            crosses = crosses + 1;
            flipflop = 1;
        end
    elseif flipflop == 1
        if data > threshold
            crosses = crosses + 1;
            flipflop = -1;
        end
    elseif flipflop == -1
        if data <  -threshold
            crosses = crosses + 1;
            flipflop = 1;
        end
    end
    
    if data > maxamp
        maxamp = data;
    end
    if data < minamp
        minamp = data;
    end
end

avg
crosses/2
minamp
maxamp
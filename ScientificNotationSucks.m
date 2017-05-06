function print = ScientificNotationSucks(Value,SigFigs)
   
if(SigFigs == 1)
    print = num2str(Value,'%.0e');
elseif(SigFigs == 2)
    print = num2str(Value,'%.1e');
elseif(SigFigs == 3)
    print = num2str(Value,'%.2e');
elseif(SigFigs == 4)
    print = num2str(Value,'%.3e');
elseif(SigFigs == 5)
    print = num2str(Value,'%.4e');
elseif(SigFigs == 6)
    print = num2str(Value,'%.5e');
elseif(SigFigs == 7)
    print = num2str(Value,'%.6e');
elseif(SigFigs == 8)
    print = num2str(Value,'%.7e');
elseif(SigFigs == 9)
    print = num2str(Value,'%.8e');
else
    print = num2str(Value,'%.9e');
end

end
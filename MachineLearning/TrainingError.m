function eps = TrainingError(Classify,Labels)

    corr = Classify == Labels;
    eps = mean(corr);
end
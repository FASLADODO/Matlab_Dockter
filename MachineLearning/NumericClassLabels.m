function [numericlabels, classes] = NumericClassLabels(labels)
%take class labels given a cell arrray and convert them to numeric labels
    classes = unique(labels);
    labs = 1:length(classes);
    labsc = strread(num2str(labs),'%s');
    numericlabels = double(nominal(labels,labsc));
end
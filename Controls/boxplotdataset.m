csatsdata = dataset('XLSFile','CSATSData.xlsx')

dataplot=[csatsdata.Warm,csatsdata.Cold,csatsdata.Expert]

figure;
title 'Warm vs. Cold vs. Expert' %%Not working

boxplot(dataplot, 'labels', {'Warm','Cold','Expert'})
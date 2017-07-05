figure

nbins = 11;

xvals = [20:3:40];

figure(1)

hist(grades,nbins)
title('Final Grades ME 5286 (mean 83)')
xlabel('Grade (out of 100)')
ylabel('Frequency')
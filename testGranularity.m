scale = 10;
max = 255;

ratioin = (scale/max)
ratioout = 1/ratioin;

sample = 1:255;

gran = sample*ratioin;

granround = round(gran);

granfull = granround*ratioout;
granfullround = round(granfull);

figure
plot(granfullround)


%running

time = [41 48]; %minutes,seconds
miles = 5.5;

timedec = time(1) + time(2)/60

minmiledec = timedec / miles;

minmile = floor(minmiledec) + (minmiledec - floor(minmiledec) )*0.6;


str = sprintf('Avg minutes per mile: %f',minmile );
disp(str)


eachmile = [1:ceil(miles)] * minmiledec;
allmiles = floor(eachmile) + (eachmile - floor(eachmile) )*0.6;

disp('Time at each mile:')
allmiles
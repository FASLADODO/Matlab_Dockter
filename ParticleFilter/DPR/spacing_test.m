origin = 11;
max = 25;
points = 13;

lin = linspace(origin,max,points);

spacing = lin(2)-lin(1)

val = 20;

idx = round( (val - origin)/spacing ) + 1

if(idx < 1)
    idx = 1;
end
if(idx > points)
   idx = points 
end

ding.poo = 1;
ding.boo = 45;
ding.party = 56;
big{1}=ding;
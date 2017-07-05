function Limit = PlotBoundingCube(Data,sigma,color,ploton)

% Data = 3 column data matrix
%sigma = percent to scale 3 sigma bounds by
%color = box color ([0 1 0])

if(nargin == 3)
   ploton = 1;
end

[NN,SS] = size(Data);
if(SS > 3)
   warning('Data matrix has more than 3 columns, ignoring extra columns')
end
if(SS < 3)
   warning('Data matrix has too few columns, NO PLOTS')
   ploton = 0;
end

%get the mean and stddev of the data
Range = DataBounds(Data);
Origin = mean(Range);
% Origin = mean(Data);
Scale = (max(Range) - min(Range))*sigma; %std(Data)*3*sigma; %3 sigma bounds

%Shift the origin
Origin = Origin - Scale./2;

%get limits of cube (works in ND
Limit = [Origin; Origin + Scale];

%NOW DO PLOTTING STUFF
if(ploton)
    % Define the vertexes of the unit cubic (hard coded for 3 dimensions)
    ver = [1 1 0;
        0 1 0;
        0 1 1;
        1 1 1;
        0 0 1;
        1 0 1;
        1 0 0;
        0 0 0];

    %  Define the faces of the unit cubic (hard coded for 3 dimensions)
    fac = [1 2 3 4;
        4 3 5 6;
        6 7 8 5;
        1 2 8 7;
        6 7 1 4;
        2 3 5 8];

    %construct the cube
    cube = [ver(:,1)*Scale(1)+Origin(1),ver(:,2)*Scale(2)+Origin(2),ver(:,3)*Scale(3)+Origin(3)];

    %plot it yo
    patch('Faces',fac,'Vertices',cube,'EdgeColor',color,'FaceColor','none','LineWidth',3);
end
end
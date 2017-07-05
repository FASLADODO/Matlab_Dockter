%Make some dumb data

nn=100;

A = mvnrnd([3,4],[1,0;0,1],nn);
B = mvnrnd([6,7],[1,0;0,1],nn);

Data = [A;B];
Labels = [ones(nn,1)*1; ones(nn,1)*-1 ];

figure
gscatter(Data(:,1),Data(:,2),Labels)

%Add in an offset term
Data = [Data,ones(length(Labels),1)];


%% Train the perceptron

[NN,SS] = size(Data);

theta = randn(SS,1)

cycles = 1000;
training_thresh = 0.0001;

%update according to equation 4
for ii = 1:cycles
    for jj = 1:NN
        %check current mapping
        f_x_theta = sign(Data*theta);
        %Update only if the sign is wrong
        if(f_x_theta(jj) ~= Labels(jj))
            theta = theta + Labels(jj).*Data(jj,:)';
        end
    end
    
    %check training error according to equation 3
    E = (1/NN)*sum(ones(NN,1) - (sign(Data*theta) == Labels) )
    if(E < training_thresh)
       break; 
    end
end

theta

slope = -theta(1)/theta(2)
offset  = -theta(3)/theta(2)

%% Visualize boundary

nv = 20;
xbound = linspace(0,max(Data(:,1)),nv)';
ybound = [xbound,ones(nv,1)]*[slope;offset];

figure
gscatter(Data(:,1),Data(:,2),Labels)
hold on
plot(xbound,ybound,'g-')
hold off
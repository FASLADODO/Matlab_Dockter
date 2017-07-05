%Train Darrens smart light data

Artery_in = load('SmartLightTrainDark.csv');
Plaque_in = load('SmartLightTrainLight.csv');

NUM_FIBERS = 63;

%get the time stamps
D_Artery.Time = Artery_in(:,1) - Artery_in(1,1);
D_Plaque.Time = Plaque_in(:,1) - Plaque_in(1,1);

%class labels
AC = 1;
PC = 2;

%placeholders
colorz = [1,2,3] + 1; %BGR add 1 for time
All_Artery = [];
All_Plaque = [];

for i = 1:NUM_FIBERS
    idxc = colorz + (i-1)*3;
    %Each fiber individually for each class
    D_Artery.Fibers{i} = Artery_in(:,idxc);
    D_Plaque.Fibers{i} = Plaque_in(:,idxc);
    
    %Each fiber individually, combined classes
    D_All{i} = [D_Artery.Fibers{i}; D_Plaque.Fibers{i} ];
    Class_All{i} = [ones(length(D_Artery.Fibers{i}),1)*AC ; ones(length(D_Plaque.Fibers{i}),1)*PC];
    
    All_Artery = [All_Artery; D_Artery.Fibers{i}];
    All_Plaque = [All_Plaque; D_Plaque.Fibers{i}];
end

%test classifier
fibnum = 25;
linclass = fitcdiscr(D_All{fibnum},Class_All{fibnum});
C1 = linclass.Coeffs(1,2).Const
W1 = linclass.Coeffs(1,2).Linear
C2 = linclass.Coeffs(2,1).Const
W2 = linclass.Coeffs(2,1).Linear
meanclass = predict(linclass,D_All{fibnum});


%% Train all parameters for each fiber

for i = 1:NUM_FIBERS
    linclass = fitcdiscr(D_All{i},Class_All{i});
    %LDA parameters
    C1 = linclass.Coeffs(1,2).Const;
    W1 = linclass.Coeffs(1,2).Linear;
    C2 = linclass.Coeffs(2,1).Const;
    W2 = linclass.Coeffs(2,1).Linear;
    
    %[W{i}, C{i}] = LDASimple(D_All{i},Class_All{i});
    W{i} = [W1';W2'];
    C{i} = [C1';C2'];
    
    %test classify
    
    %[P,Class] = LDAonline(D_All{i},W{i}, C{i});
    meanclass = predict(linclass,D_All{i});

    classify = meanclass == Class_All{i};
    
    accuracy(i) = sum(classify)/length(classify);
end



%% Does a time plot of all fibers for each class

figure
for fiber = 1:NUM_FIBERS
    scatter3(D_Artery.Fibers{fiber}(:,1),D_Artery.Fibers{fiber}(:,2),D_Artery.Fibers{fiber}(:,3),'r.');
    hold on
    scatter3(D_Plaque.Fibers{fiber}(:,1),D_Plaque.Fibers{fiber}(:,2),D_Plaque.Fibers{fiber}(:,3),'b.');
    hold off
    xlabel('blue')
    ylabel('green')
    zlabel('red')
    stit = sprintf('fiber: %d, Accuracy: %f', fiber, accuracy(fiber) );
    title(stit)
    legend('Artery','Plaque')
    axis([0 255 0 255 0 255])
    
    pause(1)
    clf
end

%% Write that shit out to file

%Plaque params first, then artery params
%const B G R
fileID = fopen('ldaparams.rod','w');

for i = 1:NUM_FIBERS
    fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f\n',C{i}(PC),W{i}(PC,1),W{i}(PC,2),W{i}(PC,3),C{i}(AC),W{i}(AC,1),W{i}(AC,2),W{i}(AC,3));
end

fclose(fileID);







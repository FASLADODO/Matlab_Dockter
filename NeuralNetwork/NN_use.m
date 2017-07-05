%%NN use

M = load('fisheriris.mat')

species_cell = {'setosa','versicolor','virginica'};
Classes = [0, 1; 
                1, 0; 
                1, 1];
colors = {[0 1 0], [1 0 0], [0 0 1] };

for ii = 1:length(M.species)
   Training{ii}.inputs = M.meas(ii,:);
   %get vector for each species
   index = find(ismember( species_cell, M.species(ii) ));
   Training{ii}.outputs = Classes(index,:);
   Training{ii}.class = index;
   %plots for funsies
   h(index) = scatter(Training{ii}.inputs(1),Training{ii}.inputs(3),20,'MarkerEdgeColor',colors{index},'MarkerFaceColor',colors{index},'LineWidth',1);
   hold on
end
hold off
title('iris data set')
xlabel('input1')
ylabel('input3')
legend([h(1),h(2),h(3)],species_cell{1},species_cell{2},species_cell{3}  )

tic
NNstruct = NN_Train(Training);
toc

%%


indy = randi([0,150])

outputNode = NN_Online(NNstruct,Training{indy}.inputs);

outputNode{end}
Training{indy}.outputs

%%

accuracy = NN_TestClassify(Training, NNstruct, Classes)
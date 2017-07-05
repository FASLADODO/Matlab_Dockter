function NetStruct = NN_Create(layers,nodes)

    for ll = 1:layers
        for jj = 1:nodes(ll)
            NetStruct.layer{ll}.node{jj}.weights = [];
            NetStruct.layer{ll}.node{jj}.threshold = 0;
        end
    end

%     NetStruct.layer{1}.node{1}.weights = [0.24,0.12,0.5,0.4];
%     NetStruct.layer{1}.node{2}.weights = [0.1,0.9,0.23,0.19];
%     NetStruct.layer{1}.node{3}.weights = [0.99,0.3,0.3,0.2];
%     NetStruct.layer{1}.node{4}.weights = [0.45,0.6,0.89,0.23];
% 
%     NetStruct.layer{2}.node{1}.weights = [0.24,0.12,0.5,0.4];
%     NetStruct.layer{2}.node{2}.weights = [0.1,0.9,0.23,0.19];
%     NetStruct.layer{2}.node{3}.weights = [0.99,0.3,0.3,0.2];
% 
%     NetStruct.layer{3}.node{1}.weights = [0.24,0.12,0.5];
%     NetStruct.layer{3}.node{2}.weights = [0.1,0.9,0.23];
%     NetStruct.layer{3}.node{3}.weights = [0.1,0.3,0.3];
% 
%     NetStruct.layer{1}.node{1}.threshold = 1;
%     NetStruct.layer{1}.node{2}.threshold = 1;
%     NetStruct.layer{1}.node{3}.threshold = 2;
%     NetStruct.layer{1}.node{4}.threshold = 2;
% 
%     NetStruct.layer{2}.node{1}.threshold = 0.5;
%     NetStruct.layer{2}.node{2}.threshold = 1;
%     NetStruct.layer{2}.node{3}.threshold = 0.8;
% 
%     NetStruct.layer{3}.node{1}.threshold = 2;
%     NetStruct.layer{3}.node{2}.threshold = 1;
%     NetStruct.layer{3}.node{3}.threshold = 2;

end
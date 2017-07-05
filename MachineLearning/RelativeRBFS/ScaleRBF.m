function scale = ScaleRBF(Data,Probability)
%integral estimate using ND area based on bounds of data

%data size = number of steps in each dimension
nn = size(Data,1)^(1/size(Data,2));

%data bounds in each dimensions
bounds = DataBounds(Data);

%get average step size of data
stepsize = diff(bounds)/nn;

%multiply each 
scale = prod(stepsize)*sum(Probability);


end
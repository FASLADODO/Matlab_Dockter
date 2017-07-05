function [V,U,N] = SingleFrenetFrame(Data)
%computes V U and N coordinate frame for a single window of data
%gets called by FrenetFrame
%V,U,N are each 3x1 Unit Vectors

V = 0; %Todo fit all data to a single vector
U = 0; %fit concavitity to second vector
N = cross(V,U);

end
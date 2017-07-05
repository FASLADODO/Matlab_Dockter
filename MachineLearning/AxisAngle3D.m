function [angle,axis,axnorm] = AxisAngle3D(X1,X2)
%compute the angle between two vectors in 3D
%also computes the axis untersecting the two vectors
%X1 = [x,y,z], X2 = [x,y,z]
%X1 and X2 should be unit vectors
%http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/

angle = acos(dot(X1,X2));
axis = cross(X1,X2);
axnorm = norm(axis);

end
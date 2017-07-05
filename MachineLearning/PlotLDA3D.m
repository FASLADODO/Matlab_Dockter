function [] = PlotLDA3D(Mdl,D,nn)

%MDL is struct from fitcdiscr()
%D is training data matrix
%nn is interpolation points

    lims = [min(D); max(D)];

    K = Mdl.Coeffs(1,2).Const;
    L = Mdl.Coeffs(1,2).Linear;
    
    %space out some vectors
    xv = linspace(lims(1,1),lims(2,1),nn);
    yv = linspace(lims(1,2),lims(2,2),nn);
    zv = linspace(lims(1,3),lims(2,3),nn);
    [xx,yy,zz] = meshgrid(xv,yv,zv);
    
    %functions
    f = @(x,y,z) K + [x y z]*L ;
    %reshapce
    v = f(xx(:),yy(:),zz(:));
    v = reshape(v,size(xx));
    
    isosurface(xx,yy,zz,v,0.1);
end
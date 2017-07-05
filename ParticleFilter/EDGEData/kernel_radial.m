function k = kernel_radial(u,type)
%simple kernel function using radial symmetric kernel for 2D
%type =
%'uniform','triangle','epanechnikov','quartic','triweight','guassian','cosine'

%See here:
%https://en.wikipedia.org/wiki/Kernel_(statistics)

%default kernel type
if(nargin == 1)
   type = 'uniform'; 
end

% check length
if(length(u) < 2)
   error('This is for 2D, U must be 1x2') 
end

switch type
    case 'uniform'
        if( norm(u) <= 1 )
            k = 0.5;
        else
            k = 0;
        end
    case 'triangle'
        if( norm(u) <= 1 )
            k = 1 - norm(u);
        else
            k = 0;
        end    
    case 'epanechnikov'
        if( norm(u) <= 1 )
            k = (3/4)*(1 - u*u');
        else
            k = 0;
        end 
    case 'quartic'
        if( norm(u) <= 1 )
            k = (15/16)*( (1 - u*u')^2);
        else
            k = 0;
        end
    case 'triweight'
        if( norm(u) <= 1 )
            k = (35/32)*( (1 - u*u')^3);
        else
            k = 0;
        end
    case 'guassian'
        if( norm(u) <= 1 )
            k = (1/sqrt(2*pi))* exp( (-1/2)*u*u' );
        else
            k = 0;
        end
    case 'cosine'
        if( norm(u) <= 1 )
            k = (pi/4)* cos( (pi/2)*mean(u) );
        else
            k = 0;
        end
    otherwise
        if( norm(u) <= 1 )
            k = (3/4)*(1 - u*u');
        else
            k = 0;
        end 
end

end
function k = kernels(u,type)
%simple kernel function 
%type =
%'uniform','triangle','epanechnikov','quartic','triweight','guassian','cosine'

%See here:
%https://en.wikipedia.org/wiki/Kernel_(statistics)

%default kernel type
if(nargin == 1)
   type = 'uniform'; 
end

% check length
if(length(u) > 1)
   error('This is for 1D, U must be 1x1') 
end

switch type
    case 'uniform'
        if( abs(u) <= 1 )
            k = 0.5;
        else
            k = 0;
        end
    case 'triangle'
        if( abs(u) <= 1 )
            k = 1 - abs(u);
        else
            k = 0;
        end    
    case 'epanechnikov'
        if( abs(u) <= 1 )
            k = (3/4)*(1 - u^2);
        else
            k = 0;
        end 
    case 'quartic'
        if( abs(u) <= 1 )
            k = (15/16)*( (1 - u^2)^2);
        else
            k = 0;
        end
    case 'triweight'
        if( abs(u) <= 1 )
            k = (35/32)*( (1 - u^2)^3);
        else
            k = 0;
        end
    case 'guassian'
        if( abs(u) <= 1 )
            k = (1/sqrt(2*pi))* exp( (-1/2)*u^2 );
        else
            k = 0;
        end
    case 'cosine'
        if( abs(u) <= 1 )
            k = (pi/4)* cos( (pi/2)*u );
        else
            k = 0;
        end
    otherwise
        if( abs(u) <= 1 )
            k = (3/4)*(1 - u^2);
        else
            k = 0;
        end 
end

end
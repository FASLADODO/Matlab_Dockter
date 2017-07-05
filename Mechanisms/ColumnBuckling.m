function [sigma_cr,sigma_allow,P_cr,P_allow,Area] = ColumnBuckling(Length,E,S_y,End_Condition,Cross_Section,Dimensions,S_F)

%Computes the critical axial stress for column buckling using different
%cross sections, end conditions and Euler vs. Johnson
%Returns critical stress, critical load, and area


%%%%%%Arguments%%%%%%
%Length = (meters)
%E = Modulus of Elasticity (Pascals)
%S_y = Yield Strength (Pascals)
%End_Condition = 'PinPin','FixFree','FixFix' or 'FixPin'
%Cross_Section = 'Circle' or 'Rectangle' or 'Hollow'
%Dimensions = [diameter], [base,height], [D_outer,D_inner] 
%S_F = safety factor


%%%%Example use
%[s_cr, S_allow, P_cr, P_allow, A] =
% ColumnBuckling(10,207*10^6,110*10^6,'FixFree','Rectangle',[1,2],3);

format long g
    
L_e = 0;

switch End_Condition
    case 'PinPin'
        L_e = Length; 
    case 'FixFree'
        L_e = 2*Length; 
    case 'FixFix'
        L_e = 0.5*Length; 
    case 'FixPin'
        L_e = 0.7*Length; 
    otherwise 
        error('End condition not recognized')
end

p = 0;
I = 0;
Area = 0;

switch Cross_Section 
    case 'Circle'
        I = (pi/64)*Dimensions(1)^4;
        Area = (pi/4)*Dimensions(1)^2;
        p = Dimensions(1)/4; % radius_gyration
    case 'Rectangle'
        I = (Dimensions(1)*Dimensions(2)^3 )/ 12;
        Area = Dimensions(1)*Dimensions(2);
        p = Dimensions(2)/2; % radius_gyration
    case 'Hollow'
        I = (pi/64)*(Dimensions(1)^4 - Dimensions(2)^4);
        Area = (pi/4)*(Dimensions(1)^2 - Dimensions(2)^2);
        p = sqrt( (Dimensions(1)^2 + Dimensions(2)^2) / 16); % radius_gyration
    otherwise 
        error('Cross Section not recognized')
end

limit = sqrt((2*pi^2*E)/S_y);

if ( (L_e/p) > limit) %Euler
    sigma_cr = (pi^2*E*I)/(Area*(L_e^2));
    disp('Euler Case')
else %Johnson
    sigma_cr = S_y - ((S_y^2)/(2*pi^2*E))*(L_e/p);
    disp('Johnson Case')
end

sigma_allow = sigma_cr/S_F;
P_cr = sigma_cr * Area;
P_allow = P_cr/S_F;

end
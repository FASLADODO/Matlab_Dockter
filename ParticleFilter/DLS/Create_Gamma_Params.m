function [phi_new] = Create_Gamma_Params(phi, gamma)

phi_new = phi;
phi_new(2) = phi(2) + phi(2)*gamma*tand(45);
phi_new(3) = phi(3) + phi(3)*gamma*tand(45);


end
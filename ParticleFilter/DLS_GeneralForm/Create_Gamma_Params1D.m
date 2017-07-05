function [phi_new] = Create_Gamma_Params1D(phi, gamma)

phi_new = phi;
phi_new(2) = phi(2) + phi(2)*gamma*tand(45);


end
function [ W, Z ] = StandardDyadCramers( alpha2, alpha3, beta2, beta3, delta2, delta3 )
    % solve for dyad (ZA and WA) with Cramers Rule Matrices
    
    denominator = [ exp(i*beta2) - 1, exp(i*alpha2) - 1;
                    exp(i*beta3) - 1, exp(i*alpha3) - 1 ];
    numerator1 = [ delta2, exp(i*alpha2) - 1;
                    delta3, exp(i*alpha3) - 1 ];
    numerator2 = [ exp(i*beta2) - 1, delta2;
                    exp(i*beta3) - 1, delta3 ];

    % Cramers rule for dyad using the 2 vector loop equations
    W = det(numerator1) / det(denominator);
    Z = det(numerator2) / det(denominator);

end


function [P_joint, P_marg_1,P_marg_2,P_cond] = ProbabilityRules(P1,P2)
    
    P_joint = P1*P2';
    
    [N1,N2] = size(P_joint);
    
    P_marg_1 = sum(P_joint,2);
    P_marg_2 = sum(P_joint,1);
    
    P_cond = P_joint.* repmat(P1,1,N2) ./ repmat(P2,1,N1); %P(Y|X)
end
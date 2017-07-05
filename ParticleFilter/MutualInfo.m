function [ H_X, H_Y, H_XY, H_YX, H_x_comma_y, I_XY ] = MutualInfo( probs )
% Solve for conditional entropy and mutual information given a matrix of
% co-propbability distributions
% Taken From chapter 2, "Elements of Information Theory", Thomas Cover

%Example Use: 
% A = [1/8,1/16,1/32,1/32;
%     1/16,1/8,1/32,1/32;
%     1/16,1/16,1/16,1/16;
%     1/4,0,0,0];
% 
% [ H_X, H_Y, H_XY, H_YX, H_x_comma_y, I_XY ] = MutualInfo( A );

    [row,col] = size(probs);
    p_x = [];
    p_y = [];
    %Marginals
    for kk = 1:col
       margx = sum(probs(:,kk) );
       p_x = [p_x, margx];
       H_Y_i(:,kk) = probs(:,kk).*(1/margx);
    end
    for kk = 1:row
       margy = sum(probs(kk,:) );
       p_y = [p_y, margy];
       H_X_i(kk,:) = probs(kk,:).*(1/margy);
    end
    
    %Marginal Entropy
    H_X = 0;
    H_Y = 0;
    for kk = 1:col
        H_X = H_X - p_x(kk)*log2(p_x(kk));
    end
    for kk = 1:row
        H_Y = H_Y - p_y(kk)*log2(p_y(kk));
    end
    
    %Entropy
    for kk = 1:row
        temphx = 0;
        for jj = 1:col
            if(H_X_i(kk,jj) ~= 0)
                temphx = temphx - H_X_i(kk,jj)*log2(H_X_i(kk,jj));
            end
        end
        H_X_yi(kk) = temphx;
    end
    for kk = 1:col
        temphy = 0;
        for jj = 1:row
            if(H_Y_i(jj,kk) ~= 0)
                temphy = temphy - H_Y_i(jj,kk)*log2(H_Y_i(jj,kk));
            end
        end
        H_Y_xi(kk) = temphy;
    end
    
    %conditional entropy
    H_XY = 0;
    for kk = 1:col
       H_XY = H_XY + p_y(kk)*H_X_yi(kk); 
    end
    
    H_YX = 0;
    for kk = 1:row
       H_YX = H_YX + p_x(kk)*H_Y_xi(kk); 
    end
    
    %Joint Entropy
    H_x_comma_y = H_X + H_YX;
    
    %Mutual Information
    I_XY = H_X - H_XY;
    I_YX = H_Y - H_YX;
    

end
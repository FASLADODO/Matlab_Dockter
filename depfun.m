function ans = depfun(var)

    ans = matlab.codetools.requiredFilesAndProducts(var);
    
    for ii = 1:length(ans)
        disp(ans(ii));
    end
end
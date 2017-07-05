function Diff = ComputeRBFDifference(pwithin,pbetween)

    %KL Divergence
%     Diff = log10(pwithin./(pbetween)); %LOG LIKELIHOOD
    Diff = pwithin.*log10(pwithin./(pbetween)); %REGULAR KL
%     Diff = (pwithin+pbetween).*log10(pwithin./(pbetween)); %NOT ACTUALLY KL 
%     Diff = pwithin.*log10(pwithin./pbetween) + pbetween.*log10(pbetween./pwithin); %SYMMETRIC KL
end
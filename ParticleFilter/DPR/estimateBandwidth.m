function H = estimateBandwidth( pdf )

sigma = std(pdf.Mu,0,2);

H = 1.06.*sigma.*pdf.n^(-(1/6))
end
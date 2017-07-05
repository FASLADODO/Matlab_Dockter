function [class,w_sums] = SkillClassificationPDF(online,pdf_nov,pdf_int,pdf_exp)

tic

display('Computing Thresholds... ')
d_thresh_nov = norm(pdf_nov.bw);
d_thresh_int = norm(pdf_int.bw);
d_thresh_exp = norm(pdf_exp.bw);

display('Computing Nearest Neighbors... ')
[index_nov,d_nov] = dsearchn(pdf_nov.pos,online);
[index_int,d_int] = dsearchn(pdf_int.pos,online);
[index_exp,d_exp] = dsearchn(pdf_exp.pos,online);

display('Making Vectors... ')
probs_at_each = [ pdf_nov.prob(index_nov),pdf_int.prob(index_int),pdf_exp.prob(index_exp) ];
d_at_each = [d_nov, d_int, d_exp];

display('Compute CumSums... ')
offset = repmat([ min(pdf_nov.bw), min(pdf_int.bw), min(pdf_exp.bw) ], length(d_at_each), 1 ); %offset to avoid divide by zero
new_d = d_at_each + offset;
w_combos = probs_at_each ./ new_d; %taking weighted combo

%take overall sum
w_sums = sum(w_combos);

[mm,class] = max(w_sums); %get class

display('Average Time ')
averageTime = toc/length(online)  % TOC, pair 1  

end
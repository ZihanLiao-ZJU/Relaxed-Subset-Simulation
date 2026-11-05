function [conditional_probabilities,log_pf,pf] = ReSuSProduct(g,conditional_probabilities,log_normalizers)
% Original Re-SuS probability-product postprocessing, Eqs. (21)/(24).
final_probability = mean(g<=0);
conditional_probabilities(end+1) = final_probability;
log_pf = log_normalizers(end)+log(final_probability);
pf = exp(log_pf);
end

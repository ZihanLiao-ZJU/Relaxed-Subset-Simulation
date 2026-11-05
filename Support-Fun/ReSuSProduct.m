function out = ReSuSProduct(raw)
    % Native probability-product estimator from the completed sampling history.
    values = cellfun(@(y)reshape(y(4,:,:), 1, []), raw.y, 'UniformOutput', false);
    normalizers = cellfun(@(g)g.log_normalizer, raw.g(1:raw.Nite));
    factors = raw.conditional_probabilities;
    final_probability = mean(values{end} <= 0);
    factors(end + 1) = final_probability;
    log_pf = normalizers(end) + log(final_probability);
    pf = exp(log_pf);
    out.estimator = 'product';
    out.conditional_probabilities = factors;
    out.sampler = raw.sampler;
    out.profile = raw.profile;
    out.samples = reshape(raw.x, 1, []);
    out.values = reshape(values, 1, []);
    out.thresholds = cellfun(@(g)[g.g_threshold, g.logpdf_threshold], raw.g(1:raw.Nite), 'UniformOutput', false);
    out.thresholds = vertcat(out.thresholds{:});
    out.log_normalizers = reshape(normalizers, 1, []);
    out.sample_counts = reshape(raw.Nsam, 1, []);
    out.n_calls = sum(raw.Ncal);
    out.converged = raw.flg_cvg;
    out.log_pf = log_pf;
    out.pf = pf;
end

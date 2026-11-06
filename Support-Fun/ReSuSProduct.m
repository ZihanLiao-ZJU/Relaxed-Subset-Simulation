function out = ReSuSProduct(raw)
    % Native probability-product estimator from the completed sampling history.
    values = cellfun(@(y)reshape(y(4,:,:), 1, []), raw.y, 'UniformOutput', false);
    % Estimate only transitions into distributions that were actually sampled.
    normalizers = zeros(1, raw.Nite);
    factors = zeros(1, raw.Nite);
    for j = 1:raw.Nite-1
        next = raw.g{j+1};
        logpdf = reshape(raw.y{j}(3,:,:), 1, []);
        factors(j) = mean(values{j} <= next.g_threshold | logpdf <= next.logpdf_threshold);
        if factors(j) <= 0
            error('ReSuS:Support', 'Intermediate distribution has no supported samples.');
        end
        normalizers(j+1) = normalizers(j) + log(factors(j));
    end
    final_probability = mean(values{end} <= 0);
    factors(end) = final_probability;
    log_pf = normalizers(end) + log(final_probability);
    pf = exp(log_pf);
    out.estimator = 'product';
    out.conditional_probabilities = factors;
    out.sampler = raw.sampler;
    out.profile = raw.intBay.profile;
    samples = cellfun(@(x)reshape(x, size(x, 1), []), raw.x, 'UniformOutput', false);
    out.samples = reshape(samples, 1, []);
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

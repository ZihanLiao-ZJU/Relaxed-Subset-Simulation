classdef ReSuS
    % Hybrid LSF-or-low-density intermediate distribution for Re-SuS.
    properties
        lsf
        Ndim
        Nfun = 4
        Nsam = 5000
        MaxIte = 100
        Verbose = false
        Ncal = 0
        Y
        G = struct('g_threshold', Inf, 'logpdf_threshold', Inf)
        profile = 'legacy'

    end

    methods

        function obj = ReSuS(g_threshold, logpdf_threshold, profile)
            if nargin > 0
                obj.G.g_threshold = g_threshold;
            end

            if nargin > 1
                obj.G.logpdf_threshold = logpdf_threshold;
            end

            if nargin > 2
                obj.profile = profile;
            end

            if ~any(strcmp(obj.profile, {'legacy', 'paper'}))
                error('ReSuS:Profile', 'Choose legacy or paper.');
            end
        end

        function obj = EvlY(obj, u)
            values = reshape(obj.lsf(u), 1, []);
            assert(numel(values) == size(u, 2) && all(isfinite(values)), 'LSF must return one finite value per column');
            obj.Y = [logGauss(u, obj.Ndim); values];
            obj.Ncal = size(u, 2);
        end

        function value = EvlLKF(obj)
            value = log(double(obj.contains(obj.Y(2,:), obj.Y(1,:))));
        end

        function value = EvlPDF(obj)
            value = obj.Y(1,:);
        end

        function mask = contains(obj, g, logpdf)
            mask = (g <= obj.G.g_threshold) | (logpdf <= obj.G.logpdf_threshold);
        end

        function obj = update(obj, g, logpdf, p0)
            n = numel(g);
            k = round(n * p0);
            if k < 1 || k >= n
                error('ReSuS:Quantile', 'round(N*p0) must be between 1 and N-1.');
            end

            gs = sort(g);
            ps = sort(logpdf);
            if strcmp(obj.profile, 'paper')
                obj.G.g_threshold = gs(k);
                obj.G.logpdf_threshold = ps(k);
            else
                obj.G.g_threshold = (gs(k) + gs(k + 1)) / 2;
                obj.G.logpdf_threshold = (ps(k) + ps(k + 1)) / 2;
            end
        end

        function yes = converged(obj)
            yes = obj.G.g_threshold <= 0;
        end

        function [xs, gs, chain_length] = select_seeds(obj, x, g, logpdf, n)
            % Preserve the supplied source's strict seed boundary and rounding.
            if strcmp(obj.profile, 'paper')
                candidates = find(obj.contains(g, logpdf));
                nc = round(n / ceil(n / numel(candidates)));
            else
                candidates = find(g < obj.G.g_threshold | logpdf < obj.G.logpdf_threshold);
                counts = round(n ./ round(n ./ (1:n)));
                nc = counts(find(counts <= numel(candidates), 1, 'last'));
            end

            if isempty(nc)
                xs = x(:, []);
                gs = g([]);
                chain_length = 0;
                return
            end

            index = candidates(randperm(numel(candidates), nc));
            xs = x(:, index);
            gs = g(index);
            chain_length = round(n / nc);
        end

    end

    methods

        function obj = start_run(obj)
            % The first samples are drawn directly from the full normal prior.
            validateattributes(obj.Ndim, {'numeric'}, {'scalar', 'integer', 'positive'});
            peak = -obj.Ndim / 2 * log(2 * pi);
            if ~(obj.G.g_threshold == Inf || obj.G.logpdf_threshold >= peak)
                error('ReSuS:InitialSupport', 'Initial hybrid domain must be verifiably the full Gaussian prior.');
            end
            if isfield(obj.G, 'L') && ~isequal(obj.G.L, eye(obj.Ndim)) && ~isequal(obj.G.L, 1)
                error('ReSuS:Reference', 'ReSuS requires the standard-normal reference.');
            end
            if isfield(obj.G, 'log_normalizer') && ~isequal(obj.G.log_normalizer, 0)
                error('ReSuS:InitialNormalizer', 'The full-prior initial log normalizer must be zero.');
            end
            obj.G.L = eye(obj.Ndim);
            obj.G.log_normalizer = 0;
        end

        function [x, y, calls, obj] = sample_level(obj, sampler, n, x_seed, g_seed, steps, ite)
            if ~isequal(obj.G.L, eye(obj.Ndim)) && ~isequal(obj.G.L, 1)
                error('ReSuS:Reference', 'ReSuS requires the standard-normal reference.');
            end
            if ite == 1
                x = randn(obj.Ndim, n);
                g = reshape(obj.lsf(x), 1, []);
                calls = size(x, 2);
            else
                lp_seed = logGauss(x_seed, obj.Ndim);
                packed = [log(double(obj.contains(g_seed, lp_seed))); lp_seed; lp_seed; g_seed];
                sampler.BayStc = obj;
                [x, packed, calls, obj] = sampler.SamGen(x_seed, packed, steps);
                g = reshape(packed(4,:,:), 1, []);
            end

            flat_x = reshape(x, obj.Ndim, []);
            if numel(g) ~= size(flat_x, 2) || any(~isfinite(g))
                error('SeMCS:LSF', 'LSF must return one finite value per column.');
            end

            lp = logGauss(flat_x, obj.Ndim);
            y = [log(double(obj.contains(g, lp))); lp; lp; g];
            y = reshape(y, 4, size(x, 2), size(x, 3));
        end

        function [x_seed, g_seed, steps, obj] = select_run_seeds(obj, x, y, previous, n)
            if ~isequal(obj.G.L, eye(obj.Ndim)) && ~isequal(obj.G.L, 1)
                error('ReSuS:Reference', 'ReSuS requires the standard-normal reference.');
            end
            values = y(4, :);
            logpdf = y(3, :);
            probability = mean(obj.contains(values, logpdf));
            if probability <= 0
                error('SeMCS:Support', 'Intermediate distribution has no supported samples.');
            end

            obj.G.log_normalizer = previous.G.log_normalizer + log(probability);
            obj.G.transition_probability = probability;
            [x_seed, g_seed, steps] = obj.select_seeds(x, values, logpdf, n);
            if isempty(g_seed)
                error('SeMCS:Seeds', 'No admissible seeds.');
            end
        end

    end
end

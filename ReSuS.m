classdef ReSuS
    % Hybrid LSF-or-low-density intermediate distribution for Re-SuS.
    properties (Constant)
        RunOptions = struct('NumSam', 5000, 'MaxIte', 100, 'Verbose', false, ...
            'MaxIteOffset', 0, 'PresetMode', 'overlay', ...
            'StopBeforeSeeds', true, 'OutputSchema', 'probability_product')
    end

    properties
        lsf
        Ndim
        Nfun = 4
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

        function obj = start_run(obj, n, max_levels, p0, states)
            % Validate and initialize only this hybrid-domain method.
            validateattributes(obj.Ndim, {'numeric'}, {'scalar', 'integer', 'positive'});
            validateattributes(n, {'numeric'}, {'scalar', 'integer', '>=', 2});
            validateattributes(max_levels, {'numeric'}, {'scalar', 'integer', 'positive'});
            validateattributes(p0, {'numeric'}, {'scalar', '>', 0, '<', 1});
            if ~isempty(states)
                obj.G = SeMCS.prescribed_state(obj.G, states{1}, true, 'overlay');
            end
            obj.G.L = eye(obj.Ndim);
            obj.G.log_normalizer = 0;
        end

        function [x, y, calls, obj] = sample_level(obj, sampler, n, x_seed, g_seed, steps, ite)
            if ite == 1
                x = randn(obj.Ndim, n);
                g = reshape(obj.lsf(x), 1, []);
                calls = size(x, 2);
            else
                lp_seed = logGauss(x_seed, obj.Ndim);
                packed = [log(double(obj.contains(g_seed, lp_seed))); lp_seed; lp_seed; g_seed];
                sampler.BayStc = obj;
                [x, packed, calls, obj] = sampler.SamGen(x_seed, packed, steps);
                x = reshape(x, obj.Ndim, []);
                g = reshape(packed(4,:,:), 1, []);
            end

            if numel(g) ~= size(x, 2) || any(~isfinite(g))
                error('SeMCS:LSF', 'LSF must return one finite value per column.');
            end

            lp = logGauss(x, obj.Ndim);
            y = [log(double(obj.contains(g, lp))); lp; lp; g];
        end

        function [x_seed, g_seed, steps, obj, probability] = select_run_seeds(obj, x, y, previous, n)
            values = y(4, :);
            logpdf = y(3, :);
            probability = mean(obj.contains(values, logpdf));
            if probability <= 0
                error('SeMCS:Support', 'Intermediate distribution has no supported samples.');
            end

            obj.G.log_normalizer = previous.G.log_normalizer + log(probability);
            [x_seed, g_seed, steps] = obj.select_seeds(x, values, logpdf, n);
            if isempty(g_seed)
                error('SeMCS:Seeds', 'No admissible seeds.');
            end
        end

    end
end

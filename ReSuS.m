classdef ReSuS
    % Hybrid LSF-or-low-density intermediate distribution for Re-SuS.
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

        function defaults = run_defaults(~)
            defaults = struct('NumSam', 5000, 'ConPrb', 0.1, 'MaxIte', 100, 'Verbose', false);
        end

        function context = start_run(~, owner, sampler)
            validateattributes(owner.IntBay.Ndim, {'numeric'}, {'scalar', 'integer', 'positive'});
            validateattributes(owner.NumSam, {'numeric'}, {'scalar', 'integer', '>=', 2});
            validateattributes(owner.MaxIte, {'numeric'}, {'scalar', 'integer', 'positive'});
            validateattributes(owner.ConPrb, {'numeric'}, {'scalar', '>', 0, '<', 1});
            context = SeMCS.start_context(owner, sampler, owner.MaxIte, 'overlay');
            context.model.G.L = eye(context.model.Ndim);
            context.model.G.log_normalizer = 0;
            context.sampler.BayStc = context.model;
            context.conditional_probabilities = [];
        end

        function [context, record] = sample_level(~, context, level)
            if level == 1
                x = randn(context.model.Ndim, context.owner.NumSam);
                g = reshape(context.model.lsf(x), 1, []);
                calls = size(x, 2);
            else
                lp_seed = logGauss(context.xs, context.model.Ndim);
                packed = [log(double(context.model.contains(context.ys, lp_seed))); lp_seed; lp_seed; context.ys];
                context.sampler.BayStc = context.model;
                [x, packed, calls, context.model] = context.sampler.SamGen(context.xs, packed, context.steps);
                x = reshape(x, context.model.Ndim, []);
                g = reshape(packed(4,:,:), 1, []);
            end

            if numel(g) ~= size(x, 2) || any(~isfinite(g))
                error('SeMCS:LSF', 'LSF must return one finite value per column.');
            end

            lp = logGauss(x, context.model.Ndim);
            y = [log(double(context.model.contains(g, lp))); lp; lp; g];
            context.g{level, 1} = context.model.G;
            context.x{level, 1} = x;
            context.y{level, 1} = y;
            context.Nsam(level, 1) = numel(g);
            context.Nsze(level,:) = [size(x, 2), 1];
            context.Ncal(level, 1) = calls;
            record = struct('x', x, 'y', y, 'values', g, 'logpdf', lp, 'calls', calls, 'count', numel(g), 'shape', [size(x, 2), 1]);
        end

        function context = advance_level(~, context, record, level)
            context.next = context.model.update(record.values, record.logpdf, context.owner.ConPrb);
            context.g{level+1, 1} = context.next.G;
        end

        function [yes, context] = run_complete(~, context, ~, phase)
            yes = false;
            if ~strcmp(phase, 'before_seeds')
                return;
            end

            if context.next.converged()
                context.complete = true;
                yes = true;
            elseif numel(context.x) == context.limit
                yes = true;
            end
        end

        function context = select_run_seeds(~, context, record, ~)
            probability = mean(context.next.contains(record.values, record.logpdf));
            if probability <= 0
                error('SeMCS:Support', 'Intermediate distribution has no supported samples.');
            end

            context.conditional_probabilities(end + 1) = probability;
            context.next.G.log_normalizer = context.model.G.log_normalizer + log(probability);
            [context.xs, context.ys, context.steps] = context.next.select_seeds(record.x, record.values, record.logpdf, context.owner.NumSam);
            if isempty(context.ys)
                error('SeMCS:Seeds', 'No admissible seeds.');
            end

            context.model = context.next;
            context.sampler.BayStc = context.model;
        end

        function out = finish_run(~, context)
            out = SeMCS.pack_output(context);
            out.conditional_probabilities = context.conditional_probabilities;
            out.sampler = class(context.sampler);
            out.profile = context.model.profile;
        end

    end
end

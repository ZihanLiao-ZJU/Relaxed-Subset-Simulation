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

    properties (Constant)
        StopBeforeSeeds = true
    end

    properties (Dependent)
        FlgCvg
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

        function x = SamGen(obj, n)
            % Direct initial sampling is valid only for the full normal prior.
            validateattributes(obj.Ndim, {'numeric'}, {'scalar', 'integer', 'positive'});
            validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive'});
            p = 0.1;
            if isfield(obj.G, 'p'), p = obj.G.p; end
            if round(n * p) < 1 || round(n * p) >= n
                error('ReSuS:Quantile', 'round(N*p0) must be between 1 and N-1.');
            end
            peak = -obj.Ndim / 2 * log(2 * pi);
            if ~(obj.G.g_threshold == Inf || obj.G.logpdf_threshold >= peak)
                error('ReSuS:InitialSupport', 'Initial hybrid domain must be verifiably the full Gaussian prior.');
            end
            factor = 1;
            if isfield(obj.G, 'L')
                factor = obj.G.L;
            elseif isfield(obj.G, 'L_pi')
                factor = obj.G.L_pi;
            end
            if ~isequal(factor, eye(obj.Ndim)) && ~isequal(factor, 1)
                error('ReSuS:Reference', 'ReSuS requires the standard-normal reference.');
            end
            x = randn(obj.Ndim, n);
        end

        function obj = EvlY(obj, u)
            factor = 1;
            if isfield(obj.G, 'L')
                factor = obj.G.L;
            elseif isfield(obj.G, 'L_pi')
                factor = obj.G.L_pi;
            end
            if ~isequal(factor, eye(obj.Ndim)) && ~isequal(factor, 1)
                error('ReSuS:Reference', 'ReSuS requires the standard-normal reference.');
            end
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

        function [y, g, obj, N_cal] = UpdObj(obj, ~, y, p)
            values = y(4,:);
            logpdf = y(3,:);
            k = round(numel(values) * p);
            if k < 1 || k >= numel(values)
                error('ReSuS:Quantile', 'round(N*p0) must be between 1 and N-1.');
            end
            gs = sort(values);
            ps = sort(logpdf);
            if strcmp(obj.profile, 'paper')
                obj.G.g_threshold = gs(k);
                obj.G.logpdf_threshold = ps(k);
            else
                obj.G.g_threshold = (gs(k) + gs(k + 1)) / 2;
                obj.G.logpdf_threshold = (ps(k) + ps(k + 1)) / 2;
            end
            g = obj.G;
            N_cal = 0;
        end

        function [xs, ys] = SedSlt(obj, x, y, n)
            % Preserve the source's profile-specific boundary and rounding.
            x = reshape(x, obj.Ndim, []);
            y = reshape(y, obj.Nfun, []);
            if strcmp(obj.profile, 'paper')
                candidates = find(obj.contains(y(4,:), y(3,:)));
                nc = round(n / ceil(n / numel(candidates)));
            else
                candidates = find(y(4,:) < obj.G.g_threshold | y(3,:) < obj.G.logpdf_threshold);
                counts = round(n ./ round(n ./ (1:n)));
                nc = counts(find(counts <= numel(candidates), 1, 'last'));
            end
            if isempty(nc) || nc <= 0
                error('SeMCS:Seeds', 'No admissible seeds.');
            end
            index = candidates(randperm(numel(candidates), nc));
            xs = x(:, index);
            ys = y(:, index);
            % Prescribed domains use the same cached model values as adaptive ones.
            ys(1,:) = log(double(obj.contains(ys(4,:), ys(3,:))));
        end

        function yes = get.FlgCvg(obj)
            yes = obj.G.g_threshold <= 0;
        end
    end
end

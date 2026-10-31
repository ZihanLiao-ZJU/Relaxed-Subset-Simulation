classdef ReSuS
    % Hybrid LSF-or-low-density intermediate distribution for Re-SuS.
    properties
        g_threshold = Inf
        logpdf_threshold = Inf
    end
    methods
        function obj = ReSuS(g_threshold, logpdf_threshold)
            if nargin > 0, obj.g_threshold = g_threshold; end
            if nargin > 1, obj.logpdf_threshold = logpdf_threshold; end
        end
        function mask = contains(obj, g, logpdf)
            mask = (g <= obj.g_threshold) | (logpdf <= obj.logpdf_threshold);
        end
        function obj = update(obj, g, logpdf, p0)
            n = numel(g); k = round(n*p0);
            if k < 1 || k >= n
                error('ReSuS:Quantile','round(N*p0) must be between 1 and N-1.');
            end
            gs = sort(g); ps = sort(logpdf);
            obj.g_threshold = (gs(k)+gs(k+1))/2;
            obj.logpdf_threshold = (ps(k)+ps(k+1))/2;
        end
        function yes = converged(obj)
            yes = obj.g_threshold <= 0;
        end
        function [xs, gs, chain_length] = select_seeds(obj, x, g, logpdf, n)
            % Preserve the supplied source's strict seed boundary and rounding.
            candidates = find(g < obj.g_threshold | logpdf < obj.logpdf_threshold);
            counts = round(n ./ round(n ./ (1:n)));
            nc = counts(find(counts <= numel(candidates),1,'last'));
            if isempty(nc)
                xs = x(:,[]); gs = g([]); chain_length = 0; return
            end
            index = candidates(randperm(numel(candidates),nc));
            xs = x(:,index); gs = g(index); chain_length = round(n/nc);
        end
    end
end

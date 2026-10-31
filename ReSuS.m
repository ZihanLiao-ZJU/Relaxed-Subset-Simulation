classdef ReSuS
    % Hybrid LSF-or-low-density intermediate distribution for Re-SuS.
    properties
        lsf
        Ndim
        Nfun = 4
        Ncal = 0
        Y
        G = struct()
        profile = 'legacy'
        g_threshold = Inf
        logpdf_threshold = Inf
    end
    methods
        function obj = ReSuS(g_threshold, logpdf_threshold, profile)
            if nargin > 0, obj.g_threshold = g_threshold; end
            if nargin > 1, obj.logpdf_threshold = logpdf_threshold; end
            if nargin > 2, obj.profile=profile; end
            if ~any(strcmp(obj.profile,{'legacy','paper'})), error('ReSuS:Profile','Choose legacy or paper.'); end
        end
        function obj = EvlY(obj,u)
            values=reshape(obj.lsf(u),1,[]);
            assert(numel(values)==size(u,2) && all(isfinite(values)), 'LSF must return one finite value per column');
            obj.Y=[-.5*sum(u.^2,1)-obj.Ndim/2*log(2*pi);values];
            obj.Ncal=size(u,2);
        end
        function value = EvlLKF(obj)
            value=log(double(obj.contains(obj.Y(2,:),obj.Y(1,:))));
        end
        function value = EvlPDF(obj)
            value=obj.Y(1,:);
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
            if strcmp(obj.profile,'paper')
                obj.g_threshold=gs(k); obj.logpdf_threshold=ps(k);
            else
                obj.g_threshold = (gs(k)+gs(k+1))/2;
                obj.logpdf_threshold = (ps(k)+ps(k+1))/2;
            end
        end
        function yes = converged(obj)
            yes = obj.g_threshold <= 0;
        end
        function [xs, gs, chain_length] = select_seeds(obj, x, g, logpdf, n)
            % Preserve the supplied source's strict seed boundary and rounding.
            if strcmp(obj.profile,'paper')
                candidates=find(obj.contains(g,logpdf));
                nc=round(n/ceil(n/numel(candidates)));
            else
                candidates = find(g < obj.g_threshold | logpdf < obj.logpdf_threshold);
                counts = round(n ./ round(n ./ (1:n)));
                nc = counts(find(counts <= numel(candidates),1,'last'));
            end
            if isempty(nc)
                xs = x(:,[]); gs = g([]); chain_length = 0; return
            end
            index = candidates(randperm(numel(candidates),nc));
            xs = x(:,index); gs = g(index); chain_length = round(n/nc);
        end
    end
end

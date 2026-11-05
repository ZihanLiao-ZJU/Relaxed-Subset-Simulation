classdef SeMCS
    % Sequential Monte Carlo orchestration with Re-SuS probability products (Eq. 24).
    properties
        lsf
        dimension
        n_samples = 5000
        p0 = 0.1
        max_levels = 100
        sampler = ES()
        profile = 'legacy'
    end
    methods
        function obj = SeMCS(lsf,dimension,n_samples,p0,max_levels,sampler,profile)
            obj.lsf = lsf; obj.dimension = dimension;
            if nargin>2, obj.n_samples = n_samples; end
            if nargin>3, obj.p0 = p0; end
            if nargin>4, obj.max_levels = max_levels; end
            if nargin>5, obj.sampler = sampler; end
            if nargin>6, obj.profile = profile; end
        end
        function out = run(obj)
            validateattributes(obj.dimension,{'numeric'},{'scalar','integer','positive'});
            validateattributes(obj.n_samples,{'numeric'},{'scalar','integer','>=',2});
            validateattributes(obj.max_levels,{'numeric'},{'scalar','integer','positive'});
            validateattributes(obj.p0,{'numeric'},{'scalar','>',0,'<',1});
            distribution = ReSuS(Inf,Inf,obj.profile);
            distribution.lsf=obj.lsf;distribution.Ndim=obj.dimension;distribution.G.L=eye(obj.dimension);
            sampler_history={};
            samples = {}; values = {}; conditional_probabilities = [];
            distributions = {}; log_normalizers = []; counts = [];
            n_calls = 0; converged = false;
            next_log_normalizer = 0;
            for level = 1:obj.max_levels
                if level == 1
                    x = randn(obj.dimension,obj.n_samples);
                    g = reshape(obj.lsf(x),1,[]); calls = size(x,2);
                else
                    lp_seed=logGauss(xseed,obj.dimension);
                    packed=[log(double(distribution.contains(gseed,lp_seed)));lp_seed;lp_seed;gseed];
                    obj.sampler.BayStc=distribution;
                    previous_state=distribution.G;
                    [x,packed,calls,distribution]=obj.sampler.SamGen(xseed,packed,chain_length);
                    sampler_history{end+1}=struct('initial',previous_state,'final',distribution.G);
                    x=reshape(x,obj.dimension,[]);g=reshape(packed(4,:,:),1,[]);
                end
                if numel(g)~=size(x,2) || any(~isfinite(g))
                    error('SeMCS:LSF','LSF must return one finite value per column.');
                end
                n_calls = n_calls+calls;
                logpdf = logGauss(x,obj.dimension);
                samples{level} = x; values{level} = g;
                distributions{level} = distribution;
                log_normalizers(level) = next_log_normalizer;
                counts(level) = numel(g);
                next = distribution.update(g,logpdf,obj.p0);
                if next.converged(), converged = true; break; end
                if level == obj.max_levels, break; end
                % Eq. (23): use this level's actual sample count, not target N.
                probability = mean(next.contains(g,logpdf));
                if probability <= 0, error('SeMCS:Support','Intermediate distribution has no supported samples.'); end
                conditional_probabilities(end+1) = probability;
                next_log_normalizer = log_normalizers(level)+log(probability);
                [xseed,gseed,chain_length] = next.select_seeds(x,g,logpdf,obj.n_samples);
                if isempty(gseed), error('SeMCS:Seeds','No admissible seeds.'); end
                distribution = next;
            end
            % Eq. (21)/(24): final true-failure factor, no extra hybrid factor.
            [conditional_probabilities,log_pf,pf] = ReSuSProduct(g,conditional_probabilities,log_normalizers);
            out.estimator = 'product'; out.conditional_probabilities = conditional_probabilities;
            out.sampler=class(obj.sampler); out.profile=obj.profile;
            out.sampler_history=sampler_history;
            out.samples = samples; out.values = values;
            out.thresholds = cellfun(@(d) [d.g_threshold,d.logpdf_threshold],distributions,'UniformOutput',false);
            out.thresholds = vertcat(out.thresholds{:});
            out.log_normalizers = log_normalizers; out.sample_counts = counts;
            out.n_calls = n_calls; out.converged = converged;
            out.log_pf = log_pf; out.pf = pf;
        end
    end
end

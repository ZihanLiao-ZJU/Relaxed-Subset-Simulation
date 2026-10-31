classdef SeMCS
    % Sequential Monte Carlo orchestration with historical balance-heuristic MIS.
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
            samples = {}; values = {}; logpdfs = {};
            distributions = {}; log_normalizers = []; counts = [];
            n_calls = 0; converged = false;
            next_log_normalizer = 0;
            for level = 1:obj.max_levels
                if level == 1
                    x = randn(obj.dimension,obj.n_samples);
                    g = reshape(obj.lsf(x),1,[]); calls = size(x,2);
                else
                    lp_seed=-.5*sum(xseed.^2,1)-obj.dimension/2*log(2*pi);
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
                logpdf = -.5*sum(x.^2,1)-obj.dimension/2*log(2*pi);
                samples{level} = x; values{level} = g; logpdfs{level} = logpdf;
                distributions{level} = distribution;
                log_normalizers(level) = next_log_normalizer;
                counts(level) = numel(g);
                next = distribution.update(g,logpdf,obj.p0);
                all_g = [values{:}]; all_logpdf = [logpdfs{:}];
                denominator = SeMCS.log_denominator(distributions,log_normalizers,counts,all_g,all_logpdf);
                if next.converged(), converged = true; break; end
                if level == obj.max_levels, break; end
                mask = next.contains(all_g,all_logpdf);
                if ~any(mask), error('SeMCS:Support','Intermediate distribution has no supported samples.'); end
                next_log_normalizer = SeMCS.logsumexp(-denominator(mask));
                [xseed,gseed,chain_length] = next.select_seeds(x,g,logpdf,obj.n_samples);
                if isempty(gseed), error('SeMCS:Seeds','No admissible seeds.'); end
                distribution = next;
            end
            log_pf = SeMCS.logsumexp(-denominator(all_g<=0));
            out.sampler=class(obj.sampler); out.profile=obj.profile;
            out.sampler_history=sampler_history;
            out.samples = samples; out.values = values;
            out.thresholds = cellfun(@(d) [d.g_threshold,d.logpdf_threshold],distributions,'UniformOutput',false);
            out.thresholds = vertcat(out.thresholds{:});
            out.log_normalizers = log_normalizers; out.sample_counts = counts;
            out.n_calls = n_calls; out.converged = converged;
            out.log_pf = log_pf; out.pf = exp(log_pf);
        end
    end
    methods(Static)
        function result = logsumexp(values)
            if isempty(values), result = -Inf; return; end
            peak = max(values);
            if isinf(peak), result = peak; else, result = peak+log(sum(exp(values-peak))); end
        end
        function denominator = log_denominator(distributions,logz,counts,g,logpdf)
            denominator = -Inf(size(g));
            for j=1:numel(distributions)
                mask = distributions{j}.contains(g,logpdf);
                term = log(counts(j))-logz(j);
                previous = denominator(mask); peak = max(previous,term);
                denominator(mask) = peak+log(exp(previous-peak)+exp(term-peak));
            end
        end
    end
end

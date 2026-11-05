classdef ReSuS
    % Hybrid LSF-or-low-density intermediate distribution for Re-SuS.
    properties
        lsf
        Ndim
        Nfun = 4
        Ncal = 0
        Y
        G = struct('g_threshold',Inf,'logpdf_threshold',Inf)
        profile = 'legacy'

    end
    methods
        function obj = ReSuS(g_threshold, logpdf_threshold, profile)
            if nargin > 0, obj.G.g_threshold = g_threshold; end
            if nargin > 1, obj.G.logpdf_threshold = logpdf_threshold; end
            if nargin > 2, obj.profile=profile; end
            if ~any(strcmp(obj.profile,{'legacy','paper'})), error('ReSuS:Profile','Choose legacy or paper.'); end
        end
        function obj = EvlY(obj,u)
            values=reshape(obj.lsf(u),1,[]);
            assert(numel(values)==size(u,2) && all(isfinite(values)), 'LSF must return one finite value per column');
            obj.Y=[logGauss(u,obj.Ndim);values];
            obj.Ncal=size(u,2);
        end
        function value = EvlLKF(obj)
            value=log(double(obj.contains(obj.Y(2,:),obj.Y(1,:))));
        end
        function value = EvlPDF(obj)
            value=obj.Y(1,:);
        end
        function mask = contains(obj, g, logpdf)
            mask = (g <= obj.G.g_threshold) | (logpdf <= obj.G.logpdf_threshold);
        end
        function obj = update(obj, g, logpdf, p0)
            n = numel(g); k = round(n*p0);
            if k < 1 || k >= n
                error('ReSuS:Quantile','round(N*p0) must be between 1 and N-1.');
            end
            gs = sort(g); ps = sort(logpdf);
            if strcmp(obj.profile,'paper')
                obj.G.g_threshold=gs(k); obj.G.logpdf_threshold=ps(k);
            else
                obj.G.g_threshold = (gs(k)+gs(k+1))/2;
                obj.G.logpdf_threshold = (ps(k)+ps(k+1))/2;
            end
        end
        function yes = converged(obj)
            yes = obj.G.g_threshold <= 0;
        end
        function [xs, gs, chain_length] = select_seeds(obj, x, g, logpdf, n)
            % Preserve the supplied source's strict seed boundary and rounding.
            if strcmp(obj.profile,'paper')
                candidates=find(obj.contains(g,logpdf));
                nc=round(n/ceil(n/numel(candidates)));
            else
                candidates = find(g < obj.G.g_threshold | logpdf < obj.G.logpdf_threshold);
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
    methods
        function defaults=run_defaults(~)
            defaults=struct('NumSam',5000,'ConPrb',0.1,'MaxIte',100,'Verbose',false);
        end
        function c=start_run(~,owner,sampler)
            validateattributes(owner.IntBay.Ndim,{'numeric'},{'scalar','integer','positive'});
            validateattributes(owner.NumSam,{'numeric'},{'scalar','integer','>=',2});
            validateattributes(owner.MaxIte,{'numeric'},{'scalar','integer','positive'});
            validateattributes(owner.ConPrb,{'numeric'},{'scalar','>',0,'<',1});
            c=SeMCS.start_context(owner,sampler,owner.MaxIte,'overlay');
            c.model.G.L=eye(c.model.Ndim);c.model.G.log_normalizer=0;c.sampler.BayStc=c.model;
            c.conditional_probabilities=[];
        end
        function [c,r]=sample_level(~,c,level)
            if level==1
                x=randn(c.model.Ndim,c.owner.NumSam);g=reshape(c.model.lsf(x),1,[]);calls=size(x,2);
            else
                lp_seed=logGauss(c.xs,c.model.Ndim);
                packed=[log(double(c.model.contains(c.ys,lp_seed)));lp_seed;lp_seed;c.ys];
                c.sampler.BayStc=c.model;
                [x,packed,calls,c.model]=c.sampler.SamGen(c.xs,packed,c.steps);
                x=reshape(x,c.model.Ndim,[]);g=reshape(packed(4,:,:),1,[]);
            end
            if numel(g)~=size(x,2)||any(~isfinite(g)),error('SeMCS:LSF','LSF must return one finite value per column.');end
            lp=logGauss(x,c.model.Ndim);y=[log(double(c.model.contains(g,lp)));lp;lp;g];
            c.g{level,1}=c.model.G;c.x{level,1}=x;c.y{level,1}=y;
            c.Nsam(level,1)=numel(g);c.Nsze(level,:)=[size(x,2),1];c.Ncal(level,1)=calls;
            r=struct('x',x,'y',y,'values',g,'logpdf',lp,'calls',calls,'count',numel(g),'shape',[size(x,2),1]);
        end
        function c=advance_level(~,c,r,level)
            c.next=c.model.update(r.values,r.logpdf,c.owner.ConPrb);
            c.g{level+1,1}=c.next.G;
        end
        function [yes,c]=run_complete(~,c,~,phase)
            yes=false;if ~strcmp(phase,'before_seeds'),return;end
            if c.next.converged(),c.complete=true;yes=true;
            elseif numel(c.x)==c.limit,yes=true;end
        end
        function c=select_run_seeds(~,c,r,~)
            probability=mean(c.next.contains(r.values,r.logpdf));
            if probability<=0,error('SeMCS:Support','Intermediate distribution has no supported samples.');end
            c.conditional_probabilities(end+1)=probability;
            c.next.G.log_normalizer=c.model.G.log_normalizer+log(probability);
            [c.xs,c.ys,c.steps]=c.next.select_seeds(r.x,r.values,r.logpdf,c.owner.NumSam);
            if isempty(c.ys),error('SeMCS:Seeds','No admissible seeds.');end
            c.model=c.next;c.sampler.BayStc=c.model;
        end
        function out=finish_run(~,c)
            out=SeMCS.pack_output(c);
            out.conditional_probabilities=c.conditional_probabilities;
            out.sampler=class(c.sampler);out.profile=c.model.profile;
        end
    end
end

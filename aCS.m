classdef aCS
    % Public-library conditional sampler: sigma0=1, persistent adaptive lambda.
    properties
        initial_lambda = 0.6
        target_acceptance = 0.44
        n_groups = 5
        current_lambda = []
        level_history = {}
    end
    methods
        function [x,g,n_calls,obj] = sample(obj,lsf,distribution,xseed,gseed,chain_length)
            [d,nc]=size(xseed); current=xseed; values=reshape(gseed,1,[]);
            x=zeros(d,nc,chain_length); g=zeros(nc,chain_length); n_calls=0;
            normals=randn(d,nc,chain_length);
            uniforms=rand(nc,chain_length); %#ok<NASGU> source indicator-MH allocation
            sigma0=1; lambda=obj.initial_lambda;
            if ~isempty(obj.current_lambda), lambda=obj.current_lambda; end
            initial_lambda=lambda; updates=struct('group',{},'acceptance',{},'final_lambda',{});
            group_size=max(ceil(nc/obj.n_groups),1);
            for group=1:obj.n_groups
                ids=(group-1)*group_size+1:min(group*group_size,nc);
                if isempty(ids), continue; end
                sigma=min(lambda*sigma0,1); rho=sqrt(1-sigma.^2); accepted_count=0;
                for step=1:chain_length
                    proposal=rho.*current(:,ids)+sigma.*normals(:,ids,step);
                    proposed_g=reshape(lsf(proposal),1,[]);
                    if numel(proposed_g)~=numel(ids) || any(~isfinite(proposed_g))
                        error('aCS:LSF','LSF must return one finite value per column.');
                    end
                    lp=-.5*sum(proposal.^2,1)-d/2*log(2*pi);
                    accepted=distribution.contains(proposed_g,lp) | (uniforms(ids,step)'==0);
                    current(:,ids(accepted))=proposal(:,accepted); values(ids(accepted))=proposed_g(accepted);
                    x(:,ids,step)=current(:,ids); g(ids,step)=values(ids)';
                    accepted_count=accepted_count+sum(accepted); n_calls=n_calls+numel(ids);
                end
                rate=accepted_count/(numel(ids)*chain_length);
                lambda=exp(log(lambda)+(rate-obj.target_acceptance)/sqrt(group));
                updates(end+1)=struct('group',group,'acceptance',rate,'final_lambda',lambda);
            end
            obj.current_lambda=lambda;
            obj.level_history{end+1}=struct('initial_lambda',initial_lambda,'final_lambda',lambda,'groups',updates);
            x=reshape(x,d,[]); g=reshape(g,1,[]);
        end
    end
end

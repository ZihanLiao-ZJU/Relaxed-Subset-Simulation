classdef HMC
    % Rejection-based exact Gaussian HMC, with paper acceptance adaptation.
    % Initial time and batch count are explicit implementation choices.
    properties
        initial_time = pi/4
        n_groups = 5
        carry_state = true
        current_time = []
    end
    methods
        function [x,g,n_calls,obj] = sample(obj,lsf,distribution,xseed,gseed,chain_length)
            [d,nc]=size(xseed); current=xseed; values=reshape(gseed,1,[]);
            x=zeros(d,nc,chain_length); g=zeros(nc,chain_length); n_calls=0;
            momenta=randn(d,nc,chain_length); time=obj.initial_time;
            if obj.carry_state && ~isempty(obj.current_time), time=obj.current_time; end
            group_size=max(ceil(nc/obj.n_groups),1);
            for group=1:obj.n_groups
                ids=(group-1)*group_size+1:min(group*group_size,nc);
                if isempty(ids), continue; end
                accepted_count=0;
                for step=1:chain_length
                    proposal=current(:,ids)*cos(time)+momenta(:,ids,step)*sin(time);
                    proposed_g=reshape(lsf(proposal),1,[]);
                    if numel(proposed_g)~=numel(ids) || any(~isfinite(proposed_g))
                        error('HMC:LSF','LSF must return one finite value per column.');
                    end
                    lp=-.5*sum(proposal.^2,1)-d/2*log(2*pi);
                    accepted=distribution.contains(proposed_g,lp);
                    current(:,ids(accepted))=proposal(:,accepted); values(ids(accepted))=proposed_g(accepted);
                    x(:,ids,step)=current(:,ids); g(ids,step)=values(ids)';
                    accepted_count=accepted_count+sum(accepted); n_calls=n_calls+numel(ids);
                end
                rate=accepted_count/(numel(ids)*chain_length);
                if rate<0.3, time=asin(min(1,sin(time)*exp((rate-0.3)/2)));
                elseif rate>0.5, time=asin(min(1,sin(time)*exp((rate-0.5)/2))); end
            end
            if obj.carry_state, obj.current_time=time; end
            x=reshape(x,d,[]); g=reshape(g,1,[]);
        end
    end
end

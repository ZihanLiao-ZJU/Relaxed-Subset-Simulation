classdef ES
    % Supplied coordinatewise angular slice sampler in standard-normal space.
    properties
        max_tries = 20
    end
    methods
        function [x,g,n_calls] = sample(obj,lsf,distribution,xseed,gseed,chain_length)
            [d,nc] = size(xseed);
            x = zeros(d,nc,chain_length); g = zeros(nc,chain_length);
            normals = randn(d,nc,chain_length);
            mh_uniforms = rand(1,nc,chain_length); %#ok<NASGU> Match public indicator-MH random allocation.
            n_calls = 0;
            for step = 1:chain_length
                current = xseed; values = gseed;
                radius = sqrt(xseed.^2 + normals(:,:,step).^2);
                ratio = xseed ./ radius; ratio(radius==0) = 0;
                original_angle = acos(max(-1,min(1,ratio)));
                pending = 1:nc; lower = zeros(d,nc); upper = pi*ones(d,nc);
                for attempt = 1:obj.max_tries
                    angle = lower + rand(d,numel(pending)).*(upper-lower);
                    proposal = cos(angle).*radius(:,pending);
                    proposed_g = reshape(lsf(proposal),1,[]);
                    if numel(proposed_g) ~= numel(pending) || any(~isfinite(proposed_g))
                        error('ES:LSF','LSF must return one finite value per column.');
                    end
                    n_calls = n_calls + numel(pending);
                    logpdf = -.5*sum(proposal.^2,1)-d/2*log(2*pi);
                    accepted = distribution.contains(proposed_g,logpdf) | (mh_uniforms(1,pending,step)==0);
                    current(:,pending(accepted)) = proposal(:,accepted);
                    values(pending(accepted)) = proposed_g(accepted);
                    rejected = ~accepted;
                    if ~any(rejected), break; end
                    pending = pending(rejected); angle = angle(:,rejected);
                    lower = lower(:,rejected); upper = upper(:,rejected);
                    above = angle > original_angle(:,pending);
                    lower(~above) = angle(~above); upper(above) = angle(above);
                end
                x(:,:,step) = current; g(:,step) = values(:);
                % Accepted states must advance even when other chains hit the cap.
                xseed = current; gseed = values;
            end
            x = reshape(x,d,[]); g = reshape(g,1,[]);
        end
    end
end

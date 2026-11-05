function logpdf = logGauss(x,dimension)
% Standard-normal log density, preserving the original operation order.
logpdf = -.5*sum(x.^2,1)-dimension/2*log(2*pi);
end

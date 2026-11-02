function out = drive_resus(example,n_samples,seed,sampler_name,profile)
% Run one Re-SuS probability-product estimate. Example: out = drive_resus('A',1000,1,'ES').
if nargin<1, example = 'A'; end
if nargin<2, n_samples = 5000; end
if nargin<3, seed = 1; end
if nargin<4, sampler_name = 'ES'; end
if nargin<5, profile='legacy'; end
switch upper(sampler_name)
    case 'ACS', sampler=aCS();
    case 'HMC', sampler=HMC();
    case 'CW_AMH', sampler=CW_aMH();
    case 'ES', sampler=ES(); % public-library coordinatewise implementation
    otherwise, error('drive_resus:Sampler','Choose CW_aMH, ES, aCS or HMC.');
end
rng(seed,'twister');
switch upper(example)
    case 'A', lsf = @LSF_2D_PiecewiseLinear;
    case 'B', lsf = @LSF_2D_ChangingTopological;
    case 'C', lsf = @LSF_2D_LinearLog;
    otherwise, error('drive_resus:Example','Choose example A, B or C.');
end
engine = SeMCS(lsf,2,n_samples,0.1,100,sampler,profile);
out = engine.run();
fprintf('%s: Pf = %.8g, levels = %d, LSF evaluations = %d, converged = %d\n', ...
    upper(example),out.pf,numel(out.sample_counts),out.n_calls,out.converged);
end

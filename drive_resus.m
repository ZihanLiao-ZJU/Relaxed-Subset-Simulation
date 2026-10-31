function out = drive_resus(example,n_samples,seed)
% Run one Re-SuS MIS estimate. Example: out = drive_resus('A',5000,1).
if nargin<1, example = 'A'; end
if nargin<2, n_samples = 5000; end
if nargin<3, seed = 1; end
rng(seed,'twister');
switch upper(example)
    case 'A', lsf = @LSF_2D_PiecewiseLinear;
    case 'B', lsf = @LSF_2D_ChangingTopological;
    otherwise, error('drive_resus:Example','Choose example A or B.');
end
engine = SeMCS(lsf,2,n_samples);
out = engine.run();
fprintf('%s: Pf = %.8g, levels = %d, LSF evaluations = %d, converged = %d\n', ...
    upper(example),out.pf,numel(out.sample_counts),out.n_calls,out.converged);
end

%% Problem and sampling settings
example = 'A';
Nsam = 5000;                  % Constant samples per level
SamplerClass = @ES;           % CW_aMH, aCS, HMC or ES
profile = 'legacy';

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'MCMC-samplers'));
addpath(fullfile(root, 'Limit-state-functions'));
addpath(fullfile(root, 'Support-Fun'));
addpath(fullfile(root, 'Operations on Logarithmic values'));

sampler = SamplerClass();

%% Define problem and run simulation
switch upper(example)
    case'A'
        lsf =@LSF_2D_PiecewiseLinear;
    case'B'
        lsf =@LSF_2D_ChangingTopological;
    case'C'
        lsf =@LSF_2D_LinearLog;
    otherwise
        error('drive_resus:Example', 'Choose example A, B or C.');
end

model = ReSuS(Inf, Inf, profile);
model.lsf = lsf;
model.Ndim = 2;
sampler.BayStc = model;
engine = SeMCS(sampler);
engine.NumSam = Nsam;
engine.ConPrb = 0.1;
engine.MaxIte = 100;
raw = engine.RunIte();
out = ReSuSProduct(raw);
fprintf('%s: Pf = %.8g, levels = %d, LSF evaluations = %d, converged = %d\n', ...
    upper(example), out.pf, numel(out.sample_counts), out.n_calls, out.converged);

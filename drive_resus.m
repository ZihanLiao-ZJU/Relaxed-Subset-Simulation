%% Performance function
model = ReSuS();
model.lsf = @LSF_2D_PiecewiseLinear;       % Example A
% model.lsf = @LSF_2D_ChangingTopological;  % Example B
% model.lsf = @LSF_2D_LinearLog;            % Example C
model.Ndim = 2;

%% Sampler and simulation
SamplerClass = @ES;                         % CW_aMH, aCS, HMC or ES
simulation = SeMCS(SamplerClass(model));
simulation.NumSam = 5000;
raw = simulation.RunIte();

%% Failure probability
out = ReSuSProduct(raw);
fprintf('%s: Pf = %.8g, levels = %d, LSF evaluations = %d, converged = %d\n', ...
    func2str(model.lsf), out.pf, numel(out.sample_counts), out.n_calls, out.converged);

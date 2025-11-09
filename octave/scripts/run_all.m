% Top-level runner for CI / local smoke tests.
addpath(genpath('octave/lib'));
addpath(genpath('octave/pipelines'));

if ~exist('data/raw', 'dir')
  mkdir('data/raw');
end
if ~exist('plots', 'dir')
  mkdir('plots');
end
if ~exist('data/processed', 'dir')
  mkdir('data/processed');
end

% Call one pipeline as a smoke test
R = pipeline_eis('data/raw', struct());
disp(R);

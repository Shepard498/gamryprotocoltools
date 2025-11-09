function run_all()
% RUN_ALL Execute all pipelines on data/raw and export to plots/ + data/processed
root = fileparts(fileparts(mfilename('fullpath'))); % .../octave
addpath(genpath(fullfile(root))); % add whole repo path
repo = fileparts(root);
dirs.raw = fullfile(repo, 'data', 'raw');
dirs.interim = fullfile(repo, 'data', 'interim');
dirs.processed = fullfile(repo, 'data', 'processed');
dirs.plots = fullfile(repo, 'plots');
ensure_dirs(dirs);


% Default options (override per-pipeline if needed)
opts = struct();
opts.svg = true; % export SVG for each figure
opts.png = false; % also export PNG
opts.fig_visible = 'on'; % 'on' to show windows; 'off' for headless
opts.decimal = '.'; % '.' or ',' (auto-detected if 'auto')
opts.timezone = 'local'; % metadata only
opts.verbose = true;


% ---- Pipelines ----
try
pipeline_eis(dirs, opts);
catch err
warning('pipeline_eis failed: %s', err.message);
end


try
pipeline_activation(dirs, opts);
catch err
warning('pipeline_activation failed: %s', err.message);
end


try
pipeline_polarizacion(dirs, opts);
catch err
warning('pipeline_polarization failed: %s', err.message);
end


try
pipeline_ocp(dirs, opts);
catch err
warning('pipeline_ocp failed: %s', err.message);
end


if opts.verbose, fprintf('Done. Exports in %s and %s\n', dirs.processed, dirs.plots); end
end


% Allow `octave -qf --eval "run('octave/scripts/run_all.m')"`
if ~isdeployed && ~ismcc
st = dbstack();
if numel(st)==1
run_all();
end
end
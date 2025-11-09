function run_all()
  root = fileparts(fileparts(mfilename('fullpath')));
  addpath(genpath(root));
  repo = fileparts(root);
  dirs.raw       = fullfile(repo, 'data', 'raw');
  dirs.interim   = fullfile(repo, 'data', 'interim');
  dirs.processed = fullfile(repo, 'data', 'processed');
  dirs.plots     = fullfile(repo, 'plots');
  ensure_dirs(dirs);

  opts = struct();
  try, pipeline_eis(dirs, opts); catch err, warning('pipeline_eis failed: %s', err.message); end
  try, pipeline_activation(dirs, opts); catch err, warning('pipeline_activation failed: %s', err.message); end
  try, pipeline_polarization(dirs, opts); catch err, warning('pipeline_polarization failed: %s', err.message); end
  try, pipeline_ocp(dirs, opts); catch err, warning('pipeline_ocp failed: %s', err.message); end
  fprintf('Done. Exports in %s and %s\n', dirs.processed, dirs.plots);
end

if ~isdeployed && ~ismcc
  st = dbstack();
  if numel(st)==1
    run_all();
  end
end


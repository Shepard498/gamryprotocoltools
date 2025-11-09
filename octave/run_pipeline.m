function R = run_pipeline(name, folder, opts)
% RUN_PIPELINE  Dispatch to a specific Gamry processing pipeline.
% Usage:
%   R = run_pipeline('activacion', 'C:\data\run1', struct('export_filename','out.xlsx'));
%
% Pipelines supported (placeholders where noted):
%   - activacion   (implemented)
%   - polarizacion (TODO)
%   - ocp          (TODO)
%   - eis          (TODO)

  if nargin < 3, opts = struct(); end
  if nargin < 2 || isempty(folder), folder = pwd; end
  if ~exist(folder, 'dir')
    error('Folder does not exist: %s', folder);
  end

  % Ensure libs/pipelines are on path
  addpath(fullfile(pwd,'lib'));
  addpath(fullfile(pwd,'pipelines'));

  % Normalize pipeline name/aliases
  pname = lower(strtrim(name));
  switch pname
    case {'activacion','activation','act'}
      pipeline_fn = @pipeline_activacion;

    case {'polarizacion','polarization','pol'}
      pipeline_fn = @pipeline_polarizacion;  % placeholder

    case {'ocp'}
      pipeline_fn = @pipeline_ocp;  % placeholder

    case {'eis'}
      pipeline_fn = @pipeline_eis;  % placeholder

    otherwise
      error('Unknown pipeline: %s', name);
  end

  % If we got here with an actual function handle, run it
  if isa(pipeline_fn, 'function_handle')
    t0 = clock;
    fprintf('[%s] Running pipeline "%s"...\n', datestr(now,'HH:MM:SS'), func2str(pipeline_fn));
    fprintf('  folder: %s\n', folder);
    R = pipeline_fn(folder, opts);
    dt = etime(clock, t0);
    fprintf('[%s] Done in %.2f s.\n', datestr(now,'HH:MM:SS'), dt);
  else
    % not_implemented() already errored above; this is just a guard
    R = struct();
  end
end

% --------- helpers ---------
function not_implemented(tag)
  error('Pipeline "%s" is not implemented yet. Choose "activacion" or add the new pipeline in /pipelines.', tag);
end


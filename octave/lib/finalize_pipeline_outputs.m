function R = finalize_pipeline_outputs(Sheets, opts, stamp, svg_jobs, export_progress)
%FINALIZE_PIPELINE_OUTPUTS  Shared tail for pipeline exports (Excel + SVG).
%   Sheets          -> cell array of {sheet_name, headers, data}
%   opts            -> pipeline options (expects export_filename, svg_dir, ...)
%   stamp           -> timestamp suffix used in SVG filenames (informational)
%   svg_jobs        -> cell array of structs with fields:
%                      handle   : figure handle to export (required)
%                      filename : full output path (required)
%                      exporter : function handle (default @export_svg_scaled)
%                      args     : cell array of extra args for exporter
%                      when     : optional logical gate (false -> skip)
%   export_progress -> optional progress fraction for workbook export.
%
%   The helper mirrors the legacy pipeline tails: it conditionally exports
%   the workbook when opts.do_export is true (or unspecified) and handles
%   SVG export with mkdir/try-catch guards when opts.save_svg is true (or
%   unspecified). Invalid figure handles are ignored silently.

  if nargin < 4 || isempty(svg_jobs)
    svg_jobs = {};
  end
  if nargin < 5 || isempty(export_progress)
    export_progress = 0.92;
  end

  R = struct('default_filename', getfield_or(opts, 'export_filename', ''), 'sheets', {Sheets});

  export_filename = getfield_or(opts, 'export_filename', R.default_filename);
  if isempty(R.default_filename), R.default_filename = export_filename; end

  % ---- Excel workbook ----
  do_export = true;
  if isfield(opts, 'do_export')
    do_export = logical(opts.do_export);
    if isempty(do_export), do_export = false; end
  end
  if do_export
    progress('Exporting workbook...', export_progress, opts);
    if isempty(export_filename)
      warning('Workbook export skipped: export_filename not provided.');
    else
      export_workbook(R, export_filename, getfield_or(opts, 'progress_cb', []));
    end
  end

  % ---- SVG exports ----
  save_svg = true;
  if isfield(opts, 'save_svg')
    save_svg = logical(opts.save_svg);
    if isempty(save_svg), save_svg = false; end
  end
  if ~save_svg || isempty(svg_jobs)
    return;
  end

  try
    svg_dir = getfield_or(opts, 'svg_dir', []);
    if isempty(svg_dir)
      error('opts.svg_dir must be defined for SVG export.');
    end
    if ~exist(svg_dir, 'dir')
      mkdir(svg_dir);
    end

    for k = 1:numel(svg_jobs)
      job = svg_jobs{k};
      if isempty(job) || ~isstruct(job)
        continue;
      end
      if isfield(job, 'when') && ~job.when
        continue;
      end
      if ~isfield(job, 'handle') || isempty(job.handle) || ~ishandle(job.handle)
        continue;
      end
      if ~isfield(job, 'filename') || isempty(job.filename)
        continue;
      end
      exporter = @export_svg_scaled;
      if isfield(job, 'exporter') && ~isempty(job.exporter)
        exporter = job.exporter;
      end
      args = {};
      if isfield(job, 'args') && ~isempty(job.args)
        args = job.args;
      end
      exporter(job.handle, job.filename, args{:});
    end
  catch ME
    warning('SVG export failed: %s', ME.message);
  end
end

function val = getfield_or(S, field, default)
  if nargin < 3, default = []; end
  if isstruct(S) && isfield(S, field)
    val = S.(field);
  else
    val = default;
  end
end

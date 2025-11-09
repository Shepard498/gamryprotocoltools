function export_workbook(R, export_filename, progress_cb)
%EXPORT_WORKBOOK  Force Excel export: one workbook, many sheets.
%   R.sheets = { {sheetname, headers_row, data_matrix_or_cell}, ... }
%   - Uses io/xlswrite *directly* (no xlsopen handle).
%   - Writes headers at A1 and data at A2 for every sheet.
%   - On failure: throws (no CSV/ODS fallback).

  if nargin < 3, progress_cb = []; end
  pb = @(msg,frac) safe_progress(progress_cb, msg, frac);

  % ---- Validate & normalize ----
  if ~isfield(R,'sheets') || ~iscell(R.sheets)
    error('R.sheets must be a cell array: { {sheet, headers, data}, ... }');
  end

  % Make sure io is available
  try
    pkg load io;
  catch err
    error('Package "io" is required for Excel export: %s', err.message);
  end

  % Normalize sheets
  for k = 1:numel(R.sheets)
    % sheet name
    if ~ischar(R.sheets{k}{1})
      error('Sheet name at sheets{%d}{1} must be char.', k);
    end
    R.sheets{k}{1} = sanitize_name(R.sheets{k}{1});

    % headers row -> 1xN cell
    hdr = R.sheets{k}{2};
    if ~iscell(hdr), error('Headers for sheet "%s" must be a cell row.', R.sheets{k}{1}); end
    if size(hdr,1) ~= 1, hdr = hdr(:)'; end
    R.sheets{k}{2} = hdr;

    % data -> cell matrix
    dat = R.sheets{k}{3};
    if isnumeric(dat)
      R.sheets{k}{3} = num2cell(dat);
    elseif iscell(dat)
      % ok
    else
      error('Data for sheet "%s" must be numeric or cell.', R.sheets{k}{1});
    end
  end

  % ---- Write all sheets into the same workbook ----
  for k = 1:numel(R.sheets)
    sh   = R.sheets{k}{1};
    head = R.sheets{k}{2};
    dat  = R.sheets{k}{3};

    pb(sprintf('Exporting "%s"...', sh), 0.90 + 0.09*k/max(1,numel(R.sheets)));

    try
      % Write headers and data using filename directly (like your older script)
      xlswrite(export_filename, head, sh, 'A1');
      if ~isempty(dat)
        xlswrite(export_filename, dat,  sh, 'A2');
      end
    catch err
      % provide helpful diagnostics and fail
      try has_jvm = usejava('jvm'); catch, has_jvm = NaN; end
      msg = sprintf(['Excel export failed for sheet "%s": %s\n' ...
                     'Diagnostics: usejava(''jvm'')=%d, Octave=%s\n' ...
                     'Make sure "pkg load io" works and Java/POI is available.'], ...
                     sh, err.message, has_jvm, version());
      error(msg);
    end
  end

  pb('Export done (XLS/XLSX).', 1.0);
end

% ---------- helpers ----------
function safe_progress(cb, msg, frac)
  if isempty(cb), return; end
  try
    if nargin < 3 || isempty(frac), frac = 0; end
    cb(msg, max(0, min(1, frac)));
  catch
    % ignore UI errors
  end
end

function s = sanitize_name(s)
  % Excel sheet name constraints: <=31 chars, no : \ / ? * [ ]
  s = regexprep(s, '[:\\/\?\*\[\]]', '_');
  if numel(s) > 31, s = s(1:31); end
  if isempty(s), s = 'Sheet1'; end
end


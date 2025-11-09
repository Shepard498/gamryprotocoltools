function R = pipeline_ocp(folder, opts)
% PIPELINE_OCP
% Files: OCP_<index>.DTA   (e.g., OCP_1.DTA, OCP_2.DTA)
%   • Imports OCP Gamry files (T, Vf, Ach) with strict column mapping:
%       Pt(1), T(2), Vf(3), Vm(4), Ach(5), Over(6), Temp(7)
%   • Plots: V vs t (overlay) and drift (dV/dt vs t)
%   • Exports one workbook with:
%       - Meta (first sheet): selected header fields from the first file
%       - One sheet PER FILE: columns = index, t_s, V_V, I_A
%       - A 'Summary' sheet with metrics per file
%   • SVG export of plots using export_svg_scaled (lib/)

  addpath(fullfile(pwd,'lib'));

  % ---------- defaults ----------
  stamp = datestr(now, 'yyyymmdd_HHMM');
  defaults = struct( ...
    'do_export', true, ...
    'export_filename', fullfile(folder, sprintf('Gamry_OCP_%s.xlsx', stamp)), ...
    'progress_cb', [], ...
    'save_svg', true, ...
    'svg_dir', fullfile(folder, 'plots_svg') ...
  );
  defaults.import = struct('decimal_comma_to_dot', true);

  if nargin < 2 || ~isstruct(opts) || isempty(opts)
    opts = struct();
  end
  opts = merge_opts(defaults, opts);

  % ---------- discover files ----------
  progress('Scanning files...', 0.04, opts);
  pat = '^OCP_(\d+)\.DTA$';
  F = list_files(folder, pat);
  if isempty(F), error('No OCP files found in: %s', folder); end

  % Parse table
  T = struct('name',{},'full',{},'fileidx',{});
  for k = 1:numel(F)
    tok = F(k).tokens;  % {fileidx}
    T(end+1).name   = F(k).name;      %#ok<AGROW>
    T(end).full     = F(k).full;
    T(end).fileidx  = str2double(tok{1});
  end
  [~,ord] = sort([T.fileidx]);  T = T(ord);

  % ---------- meta from the FIRST file ----------
  H = read_ocp_header_meta(T(1).full, opts.import);  % struct with TIMEOUT, SAMPLETIME (and a few helpful extras)
  meta_hdr = {'TIMEOUT_s','SAMPLETIME_s','AREA_cm2','FRAMEWORKVERSION','INSTRUMENTVERSION','PSTATSERIALNO'};
  meta_row = { H.TIMEOUT, H.SAMPLETIME, H.AREA, H.FRAMEWORKVERSION, H.INSTRUMENTVERSION, H.PSTATSERIALNO };

  % ---------- accumulators for summary ----------
  Sheets = {};   % { {sheetname, headers, data_cells}, ... }
  S_rows = {};   % cell rows for Summary
  S_hdr  = {'file_index','N_points','duration_s','V_mean_all','V_std_all', ...
            'V_mean_tail','V_std_tail','dVdt_mean_tail','V_mean_initial'};

  % ---------- plots prep ----------
  hVt = figure('Name','OCP: V vs t'); hold on; grid on;
  xlabel('Time (s)'); ylabel('Voltage (V)'); title('OCP: V vs t (overlay)');

  hDrift = figure('Name','OCP: dV/dt (drift)'); hold on; grid on;
  xlabel('Time (s)'); ylabel('dV/dt (V/s)'); title('OCP: drift per file');

  % ---------- loop files ----------
  for k = 1:numel(T)
    progress(sprintf('Reading %s', T(k).name), 0.10 + 0.75*(k/numel(T)), opts);

    % STRICT OCP import: T(col2), Vf(col3), Ach(col5 if present)
    [t,V,I] = import_gamry_ocp_from_table(T(k).full, opts.import);

    t = t(:); V = V(:); I = I(:);
    tl = t - t(1);                % local time starting at 0

    % per-file sheet (index equals file_index for all rows)
    headers = {'index','t_s','V_V','I_A'};
    data = [ ...
      repmat(T(k).fileidx, numel(tl), 1), ...
      tl, V, I ...
    ];
    Sheets{end+1} = {sprintf('OCP_%d', T(k).fileidx), headers, num2cell(data)}; %#ok<AGROW>

    % ---- metrics (fixed tail = last 30%) ----
    n  = numel(tl);
    if n >= 2
      dt  = diff(tl);
      dV  = diff(V) ./ dt;                          % per-sample drift
      tc  = tl(1:end-1) + 0.5*dt;                   % midpoints for dV/dt
    else
      dV  = []; tc = [];
    end
    tail_n = max(1, floor(0.30*n));
    tail_m = false(n,1); tail_m(end-tail_n+1:end) = true;
    init_n = max(1, floor(0.05*n));
    init_m = false(n,1); init_m(1:init_n) = true;

    V_mean_all   = mean(V);
    V_std_all    = std(V);
    V_mean_tail  = mean(V(tail_m));
    V_std_tail   = std(V(tail_m));
    V_mean_init  = mean(V(init_m));
    if isempty(dV)
      dVdt_mean_t = NaN;
    else
      tail_diffs = max(1, numel(dV) - (tail_n - 1) + 1) : numel(dV);
      dVdt_mean_t = mean(dV(tail_diffs));
    end

    S_rows(end+1,:) = { T(k).fileidx, n, tl(end)-tl(1), V_mean_all, V_std_all, ...
                        V_mean_tail, V_std_tail, dVdt_mean_t, V_mean_init }; %#ok<AGROW>

    % ---- plots ----
    figure(hVt);     plot(tl, V, '.-', 'LineWidth', 1.0, 'MarkerSize', 4);
    if ~isempty(dV)
      figure(hDrift); plot(tc, dV, '.', 'MarkerSize', 6);
    end
  end

  % ---------- sheets (Meta first, then per-file, then Summary) ----------
  Sheets = [ { {'Meta', meta_hdr, meta_row} }, Sheets ];   % prepend Meta
  Sheets{end+1} = {'Summary', S_hdr, S_rows};

  % ---------- export workbook ----------
  R = struct('default_filename', opts.export_filename, 'sheets', {Sheets});
  if opts.do_export
    progress('Exporting workbook...', 0.92, opts);
    export_workbook(R, opts.export_filename, opts.progress_cb);
  end

  % ---------- export SVGs ----------
  if opts.save_svg
    try
      if ~exist(opts.svg_dir,'dir'), mkdir(opts.svg_dir); end
      export_svg_scaled(hVt,    fullfile(opts.svg_dir, sprintf('OCP_V_vs_t_%s.svg', stamp)), 'Width', 1600, 'Height', 1000);
      export_svg_scaled(hDrift, fullfile(opts.svg_dir, sprintf('OCP_dVdt_%s.svg',   stamp)), 'Width', 1600, 'Height', 1000);
    catch ME
      warning('SVG export failed: %s', ME.message);
    end
  end

  progress('Done', 1.00, opts);
end

% ===== helpers (local to this pipeline) ===================================

function H = read_ocp_header_meta(filepath, imp)
% Parse a few useful header fields from an OCP DTA file.
% Returns fields as double/strings. Missing fields -> NaN / ''.
  if nargin < 2, imp = struct(); end
  if ~isfield(imp,'decimal_comma_to_dot'), imp.decimal_comma_to_dot = true; end

  H = struct('TIMEOUT',NaN,'SAMPLETIME',NaN,'AREA',NaN, ...
             'FRAMEWORKVERSION','','INSTRUMENTVERSION','','PSTATSERIALNO','');

  fid = fopen(filepath,'r','n');
  if fid < 0, error('Cannot open %s', filepath); end
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if imp.decimal_comma_to_dot, ln = strrep(ln, ',', '.'); end
    if ~isempty(strfind(ln, 'CURVE')), break; end   % stop before the table

    % QUANTs we care about
    tokQ = regexp(ln, '^([A-Z0-9_]+)\s+QUANT\s+([-\+\d\.Ee]+)\s', 'tokens');
    if ~isempty(tokQ)
      key = tokQ{1}{1}; val = str2double(tokQ{1}{2});
      switch key
        case 'TIMEOUT',     H.TIMEOUT     = val;
        case 'SAMPLETIME',  H.SAMPLETIME  = val;
        case 'AREA',        H.AREA        = val;
      end
      continue;
    end
    % LABELs we want (versions / serial)
    tokL = regexp(ln, '^([A-Z0-9_]+)\s+LABEL\s+(.+?)\s', 'tokens');
    if ~isempty(tokL)
      key = tokL{1}{1};
      switch key
        case 'FRAMEWORKVERSION',   H.FRAMEWORKVERSION = strtrim(tokL{1}{2});
        case 'INSTRUMENTVERSION',  H.INSTRUMENTVERSION = strtrim(tokL{1}{2});
        case 'PSTATSERIALNO',      H.PSTATSERIALNO     = strtrim(tokL{1}{2});
      end
    end
  end
  fclose(fid);
end

function [t,V,I] = import_gamry_ocp_from_table(filepath, imp)
  % Strict OCP reader for Gamry DTA with CURVE TABLE:
  % Columns: Pt (1), T (2), Vf (3), Vm (4), Ach (5), Over (6), Temp (7)

  if nargin < 2, imp = struct(); end
  if ~isfield(imp,'decimal_comma_to_dot'), imp.decimal_comma_to_dot = true; end

  fid = fopen(filepath,'r','n');
  if fid < 0, error('Cannot open %s', filepath); end
  raw = {};
  while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    if imp.decimal_comma_to_dot
      line = strrep(line, ',', '.');   % decimal comma → dot
    end
    raw{end+1} = line; %#ok<AGROW>
  end
  fclose(fid);

  % Collect numeric rows from the whole file (table starts after CURVE header)
  D = [];
  for j = 1:numel(raw)
    nums = sscanf(raw{j}, '%f');
    if numel(nums) >= 3
      r = NaN(1, max(7, numel(nums)));
      r(1:numel(nums)) = nums(:)';
      D(end+1,1:numel(r)) = r; %#ok<AGROW>
    end
  end
  if isempty(D), error('No numeric rows found in %s', filepath); end

  % Force the OCP mapping:
  t = D(:,2);          % T
  V = D(:,3);          % Vf
  if size(D,2) >= 5
    I = D(:,5);        % Ach (if present)
  else
    I = zeros(size(t));
  end

  % Guard against NaNs or zero-length
  bad = ~(isfinite(t) & isfinite(V));
  if any(bad), t(bad) = []; V(bad) = []; if size(D,2)>=5, I(bad)=[]; end; end
  if isempty(t), error('No valid numeric OCP samples in %s', filepath); end
end

function progress(msg, frac, opts)
  try
    if nargin < 2, frac = []; end
    if isfield(opts,'progress_cb') && ~isempty(opts.progress_cb)
      opts.progress_cb(msg, frac);
    else
      if isempty(frac), fprintf('[OCP] %s\n', msg);
      else, fprintf('[OCP] %s  (%.0f%%)\n', msg, 100*max(0,min(1,frac)));
      end
    end
  catch
  end
end


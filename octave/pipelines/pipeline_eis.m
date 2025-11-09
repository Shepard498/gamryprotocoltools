function R = pipeline_eis(folder, opts)
% PIPELINE_EIS
% Parses Gamry EIS + stabilization runs and exports sheets + plots.
% Files handled:
%   - EISPOT_0V_<idx>.DTA
%   - EISPOT_5mV_ <current>A.DTA   (comma decimal accepted, e.g. 1,5A)
%   - Est_EIS_<current>A.DTA       (chrono stabilization before EIS)
%
% Exports (Excel, same folder by default):
%   - Meta_EIS (first sheet): one row per EIS file with header params
%   - one sheet per EIS file: index, freq_Hz, Zre_ohm, Zim_ohm, |Z|_ohm, phase_deg, I_A, V_V
%   - one sheet per Est_EIS file: index, t_s, V_V, I_A
%   - Summary sheet with Rs, Rct, etc.
%
% Plots:
%   - Nyquist per file (points only) + legend   [overlay figure]
%   - Bode (|Z| and phase) + legends
%   - Nyquist overlays grouped by current (points only) + legend
%   - Rs vs current (bar chart with text tick labels; Octave-safe)
%   - Rs from 0V checks (bar chart with text tick labels; Octave-safe)
%   - Idc vs point index for each EIS file + legend
%   - NEW: One Nyquist figure **per file** (square axes)
%
% All Nyquist plots now use square axes (equal scaling).
%
% SVG export:
%   - Controlled by opts.save_svg (default true) and opts.svg_dir
%   - NEW: Per-file Nyquist SVGs saved as EIS_Nyquist_<basename>_<stamp>.svg

  addpath(fullfile(pwd,'lib'));
  clear i j

  % ---------- defaults ----------
  stamp = datestr(now, 'yyyymmdd_HHMM');
  defaults = struct( ...
    'do_export', true, ...
    'export_filename', fullfile(folder, sprintf('Gamry_EIS_%s.xlsx', stamp)), ...
    'progress_cb', [], ...
    'save_svg', true, ...
    'svg_dir', fullfile(folder, 'plots_svg'), ...
    'nyq_square', true, ...           % NEW: force square axes for Nyquist plots
    'nyq_square_pad', 0.05, ...       % NEW: 5% padding when squaring limits
    'perfile_nyq', true ...           % NEW: create one Nyquist figure per file
  );
  defaults.import = struct('decimal_comma_to_dot', true);
  if nargin < 2 || ~isstruct(opts) || isempty(opts), opts = struct(); end
  opts = merge_opts(defaults, opts);

  progress('Scanning files...', 0.04, opts);
  pat_eis0  = '^EISPOT_0V_(\d+)\.DTA$';
  pat_eisOP = '^EISPOT_5mV_([\d\.,]+)A\.DTA$';
  pat_est   = '^Est_EIS_([\d\.,]+)A\.DTA$';

  F0  = list_files(folder, pat_eis0);   % ohmic @ 0 V
  Fop = list_files(folder, pat_eisOP);  % operating-point EIS
  Fest= list_files(folder, pat_est);    % chrono stabilizations

  if isempty(F0) && isempty(Fop) && isempty(Fest)
    error('No EIS-related files found in %s', folder);
  end

  % Parse tags
  EIS = struct('name',{},'full',{},'kind',{},'idx',[],'currA',[]);
  for k=1:numel(F0)
    EIS(end+1) = struct('name',F0(k).name,'full',F0(k).full,'kind','EIS0', ...
                        'idx',str2double(F0(k).tokens{1}),'currA',NaN); %#ok<AGROW>
  end
  for k=1:numel(Fop)
    cur = strrep(Fop(k).tokens{1},',','.');
    EIS(end+1) = struct('name',Fop(k).name,'full',Fop(k).full,'kind','EISOP', ...
                        'idx',NaN,'currA',str2double(cur)); %#ok<AGROW>
  end
  % Stabilizations map by current
  EST = struct('name',{},'full',{},'currA',[]);
  for k=1:numel(Fest)
    cur = strrep(Fest(k).tokens{1},',','.');
    EST(end+1) = struct('name',Fest(k).name,'full',Fest(k).full, ...
                        'currA',str2double(cur)); %#ok<AGROW>
  end

  % Sort deterministic: EIS0 by idx, EISOP by curr
  if ~isempty(EIS)
    is0 = strcmp({EIS.kind},'EIS0');
    [~,o0]   = sort([EIS(is0).idx]);   EIS(is0) = EIS(is0)(o0);
    iop = find(~is0);
    [~,oop] = sort([EIS(iop).currA]);  EIS(iop) = EIS(iop)(oop);
  end

  % ---------- accumulators ----------
  Sheets = {};
  SumHdr = {'file','kind','curr_A','N','f_max_Hz','f_min_Hz','Rs_HF_ohm','Rlow_ohm','Rct_ohm','f_at_minIm_Hz','mag_at_minIm_ohm','Est_V_mean','Est_I_mean'};
  SumRows = {};

  % Meta per-EIS-file
  MetaHdr = {'file','kind','curr_A','idx','VDC_V','FREQINIT_Hz','FREQFINAL_Hz','PTSPERDEC','VAC_mVrms','AREA_cm2','ZGUESS_ohm','EOC_V','THD_sel','DRIFTCOR_sel','FRAMEWORKVERSION','INSTRUMENTVERSION','PSTATSERIALNO','N_points'};
  MetaRows = {};

  % Plots setup (collect handles + labels for legends)
  hNyq  = figure('Name','EIS: Nyquist (per file)'); hold on; grid on;
  xlabel('Z'' (\Omega)'); ylabel('-Z'''' (\Omega)'); title('Nyquist (points only)');
  nyq_labels = {}; nyq_handles = [];

  hBodeMag = figure('Name','EIS: |Z| vs f'); hold on; grid on;
  xlabel('f (Hz)'); ylabel('|Z| (\Omega)'); set(gca,'XScale','log'); title('Bode: magnitude');
  bodeMag_labels = {}; bodeMag_handles = [];

  hBodePh  = figure('Name','EIS: phase vs f'); hold on; grid on;
  xlabel('f (Hz)'); ylabel('Phase (deg)'); set(gca,'XScale','log'); title('Bode: phase');
  bodePh_labels = {}; bodePh_handles = [];

  % Idc vs index
  hIdcIdx = figure('Name','EIS: Idc vs point index'); hold on; grid on;
  xlabel('Point index (#)'); ylabel('I_{dc} (A)'); title('EIS: I_{dc} vs index');
  idcidx_labels = {}; idcidx_handles = [];

  % Overlay buckets by current
  Bucket = containers.Map('KeyType','char','ValueType','any'); % key '1.5A' -> struct arrays
  hGrp = []; % Will create if needed

  % NEW: per-file Nyquist figure bookkeeping + global ranges for squaring
  PerNyqFigs  = [];           % figure handles
  PerNyqNames = {};           % basename for each per-file fig
  allZre = []; allZim = [];   % for global square limits on main/grouped Nyquist

  % ---------- process EIS files ----------
  for k=1:numel(EIS)
    progress(sprintf('Reading %s', EIS(k).name), 0.10 + 0.60*(k/max(1,numel(EIS))), opts);

    % Read header meta for this file
    H = read_eis_header_meta(EIS(k).full, opts.import);

    [f, Zre, Zim, VV, II] = import_gamry_eis(EIS(k).full, opts.import);

    % ensure column vectors
    f = f(:); Zre = Zre(:); Zim = Zim(:);
    Zmag = hypot(Zre, Zim);
    phase = atan2d(Zim, Zre);  % deg, conventional

    % ---- per-file sheet ----
    headers = {'index','freq_Hz','Zre_ohm','Zim_ohm','|Z|_ohm','phase_deg','I_A','V_V'};
    data = [ (1:numel(f))', f, Zre, Zim, Zmag, phase, fill_or_nan(II, numel(f)), fill_or_nan(VV, numel(f)) ];
    basename = strip_ext(EIS(k).name);
    Sheets{end+1} = {basename, headers, num2cell(data)}; %#ok<AGROW>

    % ---- meta row for this file ----
    MetaRows(end+1,:) = { ...
      basename, EIS(k).kind, EIS(k).currA, EIS(k).idx, ...
      H.VDC, H.FREQINIT, H.FREQFINAL, H.PTSPERDEC, H.VAC, H.AREA, H.ZGUESS, H.EOC, H.THD, H.DRIFTCOR, ...
      H.FRAMEWORKVERSION, H.INSTRUMENTVERSION, H.PSTATSERIALNO, numel(f) ...
    }; %#ok<AGROW>

    % ---- metrics ----
    N = numel(f);
    if N >= 5
      M = min(8, N);
      [~, idxHF] = sort(f, 'descend');
      hf = idxHF(1:M);
      p = polyfit(-Zim(hf), Zre(hf), 1); % fit Zre vs -Zim near HF
      Rs = polyval(p, 0);                % intercept at -Zim=0
    else
      Rs = NaN;
    end
    [~, iMinIm] = min(Zim);
    fMinIm   = f(iMinIm);
    magMinIm = Zmag(iMinIm);
    Rlow     = Zre( argmin(f) );         % lowest frequency Re(Z)
    Rct      = Rlow - Rs;

    % Try to attach stabilization means (match by current if EISOP)
    estV = NaN; estI = NaN;
    if strcmp(EIS(k).kind,'EISOP') && ~isnan(EIS(k).currA)
      keyA = EIS(k).currA;
      [estV, estI] = find_est_means(EST, keyA, opts.import);
    end

    SumRows(end+1,:) = { basename, EIS(k).kind, EIS(k).currA, N, max(f), min(f), Rs, Rlow, Rct, fMinIm, magMinIm, estV, estI }; %#ok<AGROW>

    % ---- plots ----
    % Overlay Nyquist
    figure(hNyq);
    h = plot(Zre, -Zim, '.', 'MarkerSize', 14);
    nyq_handles(end+1) = h; %#ok<AGROW>
    nyq_labels{end+1}  = basename; %#ok<AGROW>

    % Collect for global square limits
    allZre = [allZre; Zre(:)];
    allZim = [allZim; Zim(:)];

    % Bode
    figure(hBodeMag);
    hm = plot(f, Zmag, '.', 'MarkerSize', 14);
    bodeMag_handles(end+1) = hm; %#ok<AGROW>
    bodeMag_labels{end+1}  = basename; %#ok<AGROW>

    figure(hBodePh);
    hp = plot(f, phase, '.', 'MarkerSize', 14);
    bodePh_handles(end+1) = hp; %#ok<AGROW>
    bodePh_labels{end+1}  = basename; %#ok<AGROW>

    % Idc vs point index (only if Idc exists)
    if ~isempty(II) && ~all(isnan(II))
      figure(hIdcIdx);
      hi = plot(1:numel(II), II(:), '.-', 'MarkerSize', 14);
      idcidx_handles(end+1) = hi; %#ok<AGROW>
      idcidx_labels{end+1}   = basename; %#ok<AGROW>
    end

    % overlay buckets
    if strcmp(EIS(k).kind,'EISOP') && ~isnan(EIS(k).currA)
      lab = sprintf('%.3gA', EIS(k).currA);
      if ~isKey(Bucket, lab), Bucket(lab) = {}; end
      val = Bucket(lab);
      val{end+1} = struct('Zre',Zre,'Zim',Zim); %#ok<AGROW>
      Bucket(lab) = val;
    end

    % ---------- NEW: per-file Nyquist figure ----------
    if opts.perfile_nyq
      hf = figure('Name', ['Nyquist: ' basename]); hold on; grid on; box on;
      xlabel('Z'' (\Omega)'); ylabel('-Z'''' (\Omega)'); title(['Nyquist - ' basename],'Interpreter','none');
      plot(Zre, -Zim, '.', 'MarkerSize', 14);
      if opts.nyq_square
        set_nyquist_square(gca, Zre, -Zim, opts.nyq_square_pad);
      end
      PerNyqFigs(end+1)  = hf;           %#ok<AGROW>
      PerNyqNames{end+1} = basename;     %#ok<AGROW>
    end
  end

  % Attach legends
  if ~isempty(nyq_handles),     figure(hNyq);     legend(nyq_handles, nyq_labels, 'Interpreter','none'); end
  if ~isempty(bodeMag_handles), figure(hBodeMag); legend(bodeMag_handles, bodeMag_labels, 'Interpreter','none'); end
  if ~isempty(bodePh_handles),  figure(hBodePh);  legend(bodePh_handles,  bodePh_labels,  'Interpreter','none'); end
  if ~isempty(idcidx_handles),  figure(hIdcIdx);  legend(idcidx_handles,  idcidx_labels,  'Interpreter','none'); end

  % Square axes for the main Nyquist overlay
  if opts.nyq_square && ~isempty(allZre)
    figure(hNyq);
    set_nyquist_square(gca, allZre, -allZim, opts.nyq_square_pad);
  end

  % ---------- overlay by current (Nyquist groups) ----------
  if ~isempty(keys(Bucket))
    hGrp = figure('Name','Nyquist overlays by current'); hold on; grid on; box on;
    xlabel('Z'' (\Omega)'); ylabel('-Z'''' (\Omega)'); title('Nyquist grouped by current');
    labs = keys(Bucket);
    for i=1:numel(labs)
      cellarr = Bucket(labs{i});
      for j=1:numel(cellarr)
        plot(cellarr{j}.Zre, -cellarr{j}.Zim, '.-', 'MarkerSize', 14);
      end
    end
    legend(labs);
    % Square axes for grouped plot (use global ranges if available)
    if opts.nyq_square && ~isempty(allZre)
      set_nyquist_square(gca, allZre, -allZim, opts.nyq_square_pad);
    end
  end

  % ---------- process stabilization files (export + Vdc vs t plot) ----------
  if ~isempty(EST)
    hEst = figure('Name','Stabilizations: Vdc vs t'); hold on; grid on;
    xlabel('t (s)'); ylabel('V_{dc} (V)'); title('Est_EIS: Vdc vs t');
    est_labels = {}; est_handles = [];

    for k=1:numel(EST)
      [t,V,I] = import_gamry_chrono_flexible(EST(k).full, opts.import); %#ok<ASGLU>
      tl = t(:) - t(1);
      headers = {'index','t_s','V_V','I_A'};
      dat = [ (1:numel(tl))', tl, V(:), I(:) ];
      Sheets{end+1} = {strip_ext(EST(k).name), headers, num2cell(dat)}; %#ok<AGROW>
      figure(hEst);
      he = plot(tl, I(:), '.-', 'MarkerSize', 14);
      est_handles(end+1) = he; %#ok<AGROW>
      est_labels{end+1}  = strip_ext(EST(k).name); %#ok<AGROW>
    end
    if ~isempty(est_handles), legend(est_handles, est_labels, 'Interpreter','none'); end
  end

  % ---------- Rs bar charts (Octave-safe tick labels) ----------
  % 1) Rs vs current (EISOP)
  C = []; RS = []; LabC = {};
  for r=1:size(SumRows,1)
    if strcmp(SumRows{r,2},'EISOP') && ~isnan(SumRows{r,3}) && ~isnan(SumRows{r,7})
      C(end+1)   = SumRows{r,3}; %#ok<AGROW>
      RS(end+1)  = SumRows{r,7}; %#ok<AGROW>
    end
  end
  if ~isempty(C)
    [C, ord] = sort(C); RS = RS(ord);
    for i=1:numel(C)
      LabC{i} = sprintf('%.3g A', C(i)); %#ok<AGROW>
    end
    hRs = figure('Name','Rs vs current (bar)'); grid on;
    n = numel(RS);
    bar(1:n, RS);
    xlim([0.5 n+0.5]);
    set(gca,'XTick',1:n,'XTickLabel',LabC);
    ylabel('R_s (\Omega)'); title('High-frequency intercept vs current');
  end

  % 2) Rs for 0V checks (EIS0)
  RS0 = []; L0 = {};
  for r=1:size(SumRows,1)
    if strcmp(SumRows{r,2},'EIS0') && ~isnan(SumRows{r,7})
      RS0(end+1) = SumRows{r,7}; %#ok<AGROW>
      L0{end+1}  = SumRows{r,1}; %#ok<AGROW> % file base name
    end
  end
  if ~isempty(RS0)
    hRs0 = figure('Name','Rs from 0V checks (bar)'); grid on;
    n0 = numel(RS0);
    bar(1:n0, RS0);
    xlim([0.5 n0+0.5]);
    set(gca,'XTick',1:n0,'XTickLabel',L0);
    ylabel('R_s (\Omega)'); title('Ohmic resistance from EISPOT_0V');
  end

  % ---------- export ----------
  % Prepend Meta_EIS sheet
  Sheets = [ { {'Meta_EIS', MetaHdr, MetaRows} }, Sheets ];
  Sheets{end+1} = {'Summary', SumHdr, SumRows};

  svg_jobs = {
    struct('handle', hNyq, 'filename', fullfile(opts.svg_dir, sprintf('EIS_Nyquist_%s.svg', stamp)), ...
           'args', {{'Width', 1600, 'Height', 1000}}), ...
    struct('handle', hBodeMag, 'filename', fullfile(opts.svg_dir, sprintf('EIS_BodeMag_%s.svg', stamp)), ...
           'args', {{'Width', 1600, 'Height', 1000}}), ...
    struct('handle', hBodePh, 'filename', fullfile(opts.svg_dir, sprintf('EIS_BodePh_%s.svg', stamp)), ...
           'args', {{'Width', 1600, 'Height', 1000}})
  };
  svg_jobs{end+1} = struct('handle', hIdcIdx, 'filename', fullfile(opts.svg_dir, sprintf('EIS_Idc_vs_index_%s.svg', stamp)), ...
                           'args', {{'Width', 1600, 'Height', 1000}}, 'when', ~isempty(idcidx_handles));
  if exist('hEst','var')
    svg_jobs{end+1} = struct('handle', hEst, 'filename', fullfile(opts.svg_dir, sprintf('EIS_Stabilizations_Vdc_vs_t_%s.svg', stamp)), ...
                             'args', {{'Width', 1600, 'Height', 1000}}, 'when', ~isempty(EST));
  end
  if exist('hGrp','var')
    svg_jobs{end+1} = struct('handle', hGrp, 'filename', fullfile(opts.svg_dir, sprintf('EIS_Nyquist_grouped_%s.svg', stamp)), ...
                             'args', {{'Width', 1600, 'Height', 1000}});
  end
  if exist('hRs','var')
    svg_jobs{end+1} = struct('handle', hRs, 'filename', fullfile(opts.svg_dir, sprintf('EIS_Rs_vs_current_%s.svg', stamp)), ...
                             'args', {{'Width', 1600, 'Height', 1000}});
  end
  if exist('hRs0','var')
    svg_jobs{end+1} = struct('handle', hRs0, 'filename', fullfile(opts.svg_dir, sprintf('EIS_Rs_0V_%s.svg', stamp)), ...
                             'args', {{'Width', 1600, 'Height', 1000}});
  end
  if ~isempty(PerNyqFigs)
    for i=1:numel(PerNyqFigs)
      if ishandle(PerNyqFigs(i))
        svg_jobs{end+1} = struct('handle', PerNyqFigs(i), ...
          'filename', fullfile(opts.svg_dir, sprintf('EIS_Nyquist_%s_%s.svg', PerNyqNames{i}, stamp)), ...
          'args', {{'Width', 1200, 'Height', 1200}});
      end
    end
  end

  R = finalize_pipeline_outputs(Sheets, opts, stamp, svg_jobs, 0.92);

  progress('Done', 1.00, opts);
end

% ================= helpers =================

function set_nyquist_square(ax, x, y, padfrac)
% Square Nyquist axes: equal scaling and same span on both axes with padding.
  if nargin<4 || isempty(padfrac), padfrac = 0.05; end
  if isempty(x) || isempty(y) || all(~isfinite(x)) || all(~isfinite(y)), return; end
  x = x(isfinite(x)); y = y(isfinite(y));
  if isempty(x) || isempty(y), return; end
  xmin = min(x); xmax = max(x);
  ymin = min(y); ymax = max(y);
  % symmetric span around data centers, then equalize
  cx = 0.5*(xmin+xmax); cy = 0.5*(ymin+ymax);
  rx = 0.5*(xmax-xmin); ry = 0.5*(ymax-ymin);
  r  = max(rx, ry);
  if r <= 0, r = max(1, abs(cx)+abs(cy))*0.01; end
  r  = r * (1+padfrac);
  xlim(ax, [cx - r, cx + r]);
  ylim(ax, [cy - r, cy + r]);
  axis(ax, 'equal');
  box(ax, 'on');
end

function H = read_eis_header_meta(filepath, imp)
% Parse key EIS header fields from a DTA file up to ZCURVE.
% Normalizes decimal commas to dots. Missing -> NaN or ''.

  if nargin<2 || ~isstruct(imp), imp = struct(); end
  if ~isfield(imp,'decimal_comma_to_dot'), imp.decimal_comma_to_dot = true; end

  H = struct('VDC',NaN,'FREQINIT',NaN,'FREQFINAL',NaN,'PTSPERDEC',NaN,'VAC',NaN,'AREA',NaN, ...
             'ZGUESS',NaN,'EOC',NaN,'THD',NaN,'DRIFTCOR',NaN, ...
             'FRAMEWORKVERSION','','INSTRUMENTVERSION','','PSTATSERIALNO','');

  fid = fopen(filepath,'r','n');
  if fid < 0, error('Cannot open %s', filepath); end
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if imp.decimal_comma_to_dot, ln = strrep(ln, ',', '.'); end
    if ~isempty(strfind(ln, 'ZCURVE')), break; end % stop before table

    % POTEN (VDC)
    tokP = regexp(ln, '^([A-Z0-9_]+)\s+POTEN\s+([-\+\d\.Ee]+)\s', 'tokens');
    if ~isempty(tokP)
      key = tokP{1}{1}; val = str2double(tokP{1}{2});
      if strcmp(key,'VDC'), H.VDC = val; end
      continue;
    end

    % QUANTs
    tokQ = regexp(ln, '^([A-Z0-9_]+)\s+QUANT\s+([-\+\d\.Ee]+)\s', 'tokens');
    if ~isempty(tokQ)
      key = tokQ{1}{1}; val = str2double(tokQ{1}{2});
      switch key
        case 'FREQINIT',   H.FREQINIT   = val;
        case 'FREQFINAL',  H.FREQFINAL  = val;
        case 'PTSPERDEC',  H.PTSPERDEC  = val;
        case 'VAC',        H.VAC        = val;
        case 'AREA',       H.AREA       = val;
        case 'ZGUESS',     H.ZGUESS     = val;
        case 'EOC',        H.EOC        = val;
      end
      continue;
    end

    % SELECTOR values we care about
    tokS = regexp(ln, '^([A-Z0-9_]+)\s+SELECTOR\s+([-\+\d\.Ee]+)\s', 'tokens');
    if ~isempty(tokS)
      key = tokS{1}{1}; val = str2double(tokS{1}{2});
      switch key
        case 'THD',      H.THD      = val;
        case 'DRIFTCOR', H.DRIFTCOR = val;
      end
      continue;
    end

    % LABELs for traceability
    tokL = regexp(ln, '^([A-Z0-9_]+)\s+LABEL\s+(.+?)\s', 'tokens');
    if ~isempty(tokL)
      key = tokL{1}{1}; sval = strtrim(tokL{1}{2});
      switch key
        case 'FRAMEWORKVERSION',   H.FRAMEWORKVERSION  = sval;
        case 'INSTRUMENTVERSION',  H.INSTRUMENTVERSION = sval;
        case 'PSTATSERIALNO',      H.PSTATSERIALNO     = sval;
      end
      continue;
    end
  end
  fclose(fid);
end

function [f, Zre, Zim, V, I] = import_gamry_eis(filepath, imp)
% Header-aware EIS import for Gamry ZCURVE tables.

  if nargin<2 || ~isstruct(imp), imp = struct(); end
  if ~isfield(imp,'decimal_comma_to_dot'), imp.decimal_comma_to_dot = true; end

  % --- read file as lines, normalizing decimal commas ---
  fid = fopen(filepath,'r','n');
  assert(fid>=0, 'Cannot open %s', filepath);
  C = {};
  while ~feof(fid)
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if imp.decimal_comma_to_dot, ln = strrep(ln, ',', '.'); end
    C{end+1} = ln; %#ok<AGROW>
  end
  fclose(fid);

  % --- find ZCURVE block + header line ---
  iZ = find(startsWith(C, 'ZCURVE'), 1, 'first');
  assert(~isempty(iZ), 'ZCURVE block not found in %s', filepath);

  % The next two lines are the column names and the units.
  hdrLine  = strtrim(C{iZ+1});
  % unitLine = strtrim(C{iZ+2}); %#ok<NASGU>

  % Tokenize header names (tab or multiple spaces)
  toks = regexp(hdrLine, '\s+', 'split');
  if isempty(toks{1}), toks(1) = []; end   % some files have a leading tab

  % Build a name->index map
  idxMap = containers.Map('KeyType','char','ValueType','double');
  for j=1:numel(toks), idxMap(toks{j}) = j; end

  % Columns we care about (case-sensitive as in your sample)
  need = {'Freq','Zreal','Zimag'};
  for j=1:numel(need)
    assert(isKey(idxMap, need{j}), 'Column "%s" not found in %s', need{j}, filepath);
  end
  hasIdc = isKey(idxMap,'Idc');  hasVdc = isKey(idxMap,'Vdc');

  % --- read numeric rows until table ends ---
  D = [];
  for k = iZ+3 : numel(C)
    if isempty(strtrim(C{k})), break; end
    nums = sscanf(C{k}, '%f');
    if isempty(nums), break; end
    % Align to last numel(toks) values (rows may start with an index)
    if numel(nums) >= numel(toks)
      row = nums(end-numel(toks)+1:end).';
      D(end+1,1:numel(row)) = row; %#ok<AGROW>
    end
  end
  assert(~isempty(D), 'No numeric ZCURVE rows parsed in %s', filepath);

  % --- extract data ---
  f   = D(:, idxMap('Freq'));
  Zre = D(:, idxMap('Zreal'));
  Zim = D(:, idxMap('Zimag'));
  if hasIdc, I = D(:, idxMap('Idc')); else, I = NaN(size(f)); end
  if hasVdc, V = D(:, idxMap('Vdc')); else, V = NaN(size(f)); end
end

function [t,V,I] = import_gamry_chrono_flexible(filepath, imp)
% Flexible chrono import for Est_EIS_* (t,V,I expected, but robust).
  [D,~] = read_numeric_table(filepath, imp.decimal_comma_to_dot);
  nc = size(D,2);
  score = -inf(1,nc);
  for c=1:nc
    v = D(:,c);
    inc = mean(diff(v) > 0);
    span = max(v) - min(v);
    score(c) = 2*inc + span/(max(1,abs(mean(v))));
  end
  [~, tcol] = max(score);
  t = D(:,tcol);

  rem = setdiff(1:nc, tcol);
  vals = abs(median(D(:,rem)));
  varv = var(D(:,rem));
  [~,vIdx] = max(varv ./ (1+abs(vals-2)));
  V = D(:, rem(vIdx));

  rem2 = setdiff(rem, rem(vIdx));
  if isempty(rem2)
    I = zeros(size(t));
  else
    [~,iIdx] = max(var(D(:,rem2)));
    I = D(:, rem2(iIdx));
  end
end

function [Vmean, Imean] = find_est_means(EST, currA, imp)
  Vmean = NaN; Imean = NaN;
  if isempty(EST), return; end
  diffs = arrayfun(@(e) abs(e.currA - currA), EST);
  [m, ix] = min(diffs);
  if isempty(m) || m > 1e-3, return; end
  [t,V,I] = import_gamry_chrono_flexible(EST(ix).full, imp); %#ok<ASGLU>
  n = numel(V);
  tail = max(1, floor(0.30*n)):n;
  Vmean = mean(V(tail));
  Imean = mean(I(tail));
end

function [D, hadComma] = read_numeric_table(filepath, fixComma)
  fid = fopen(filepath,'r','n');
  if fid < 0, error('Cannot open %s', filepath); end
  raw = {};
  while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    if fixComma, line = strrep(line, ',', '.'); end
    raw{end+1} = line; %#ok<AGROW>
  end
  fclose(fid);
  hadComma = fixComma;

  D = [];
  for j=1:numel(raw)
    nums = sscanf(raw{j}, '%f');
    if ~isempty(nums)
      D(end+1,1:numel(nums)) = nums(:)'; %#ok<AGROW>
    end
  end
  if isempty(D), error('No numeric table found in %s', filepath); end
end

function v = fill_or_nan(x, n)
  if nargin<2, n = numel(x); end
  if isempty(x), v = nan(n,1);
  elseif numel(x)==n, v = x(:);
  else v = nan(n,1);
  end
end

function b = strip_ext(fn)
  [~,b,~] = fileparts(fn);
end

function idx = argmin(v)
  [~,idx] = min(v);
end


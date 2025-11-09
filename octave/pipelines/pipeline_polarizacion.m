function R = pipeline_polarizacion(folder, opts)
% PIPELINE_POLARIZACION — robust step split + last/alpha picks, duplicate-step warnings
% Files: Curva_Polarizacion_(Asc|Dsc)_(dwell)_(proto)_#<fileidx>.DTA
%
% Each file contains TWO current steps.
% We:
%   1) Import t, V, I and read header QUANTs (SAMPLETIME, TSTEP1, TSTEP2, ISTEP1, ISTEP2...).
%   2) Detect the step boundary purely from I(t): largest absolute jump in a
%      lightly smoothed current trace (movmedian), excluding a small edge guard.
%   3) Define steps by index:
%        step1 = 1:k         (alpha index between first/last)
%        step2 = k+1:N       (alpha index between first/last)
%   4) Exports:
%        - All_Asc / All_Dsc    : every sample (index, fileidx, step, t_local, V, I, filename)
%        - IV_last_Asc / _Dsc   : exactly two rows per file (last sample of each step)
%        - IV_alpha_Asc / _Dsc  : alpha-picked sample of each step (alpha in [0,1])
%        - Meta                 : representative header values (for convenience)
%        - Per_file_QA          : counts and boundary diagnostics
%   5) Plots:
%        - V vs I (all samples, Asc/Dsc)
%        - IV (only the two LAST points per file, Asc/Dsc)
%        - NEW: V vs t (overlay Asc+Dsc) — exported to SVG
%        - NEW: I vs t (overlay Asc+Dsc) — exported to SVG
%   6) Warnings:
%        - Detect duplicate step currents across Dsc files (middle points duplicated)
%          and issue a warning with the repeated setpoint values.

  addpath(fullfile(pwd,'lib'));

  % ---------- defaults ----------
  stamp = datestr(now, 'yyyymmdd_HHMMSS');  % seconds to ensure unique filenames
  defaults = struct( ...
    'do_export', true, ...
    'export_filename', fullfile(folder, sprintf('Gamry_Polarizacion_%s.xlsx', stamp)), ...
    'progress_cb', [], ...
    'save_svg', true, ...
    'svg_dir', fullfile(folder, 'plots_svg'), ...
    'alpha', 1.0, ...                              % 0 = first, 1 = last (discrete at edges)
    'alpha_round', 'nearest', ...                  % 'nearest'|'floor'|'ceil' (for interior alphas)
    'diff_smooth_pts', 5, ...                      % NEW: odd window for smoothing V(I) before dV/dI (1 = off)
    'dup_tol', 1e-6 ...                            % Tolerance when comparing current duplicates
  );
  defaults.import = struct('header_guess_lines',160,'min_numeric_cols',5,'decimal_comma_to_dot',true);
  defaults.detect = struct('medwin',5,'edge_guard_frac',0.06);  % ~5 samples, ignore ~6% at each end

  if nargin < 2, opts = struct(); end
  opts = merge_opts(defaults, opts);
  % clamp alpha
  if ~isfield(opts,'alpha') || isempty(opts.alpha), opts.alpha = 1.0; end
  opts.alpha = max(0,min(1,opts.alpha));

  % ---------- discover files ----------
  progress('Scanning files...', 0.04, opts);
  pat = '^Curva_Polarizacion_(Asc|Dsc)_(.+?)_(.+?)_#(\d+)\.DTA$';
  F = list_files(folder, pat);
  if isempty(F), error('No polarization files found in: %s', folder); end

  % parse into struct table
  T = struct('name',{},'full',{},'dir',{},'dwell',{},'proto',{},'fileidx',{});
  for k = 1:numel(F)
    tok = F(k).tokens; % {dir,dwell,proto,fileidx}
    T(end+1).name   = F(k).name;      %#ok<AGROW>
    T(end).full     = F(k).full;
    T(end).dir      = tok{1};
    T(end).dwell    = tok{2};
    T(end).proto    = tok{3};
    T(end).fileidx  = str2double(tok{4});
  end

  % split by direction and sort by fileidx
  A = T(strcmp({T.dir},'Asc'));
  D = T(strcmp({T.dir},'Dsc'));
  [~,ia] = sort([A.fileidx]); A = A(ia);
  [~,id] = sort([D.fileidx]); D = D(id);

  % ---------- accumulators ----------
  Asc_all = []; Asc_names = {};
  Dsc_all = []; Dsc_names = {};
  % columns: [index, fileidx, step_in_file, t_local_s, V_V, I_A]

  Asc_last = []; Asc_last_names = {};
  Dsc_last = []; Dsc_last_names = {};
  % columns: [fileidx, step_in_file, t_pick_s, V_pick, I_pick]

  Asc_alpha = []; Asc_alpha_names = {};
Dsc_alpha = []; Dsc_alpha_names = {};

  % --- Global-time accumulators for plotting (per direction) ---
  asc_offset = 0; dsc_offset = 0;  % seconds
  Asc_tglob = []; Asc_Vglob = []; Asc_Iglob = [];
  Dsc_tglob = []; Dsc_Vglob = []; Dsc_Iglob = [];

  meta_rep = [];        % representative header values from first file seen
  QA_rows = {};         % diagnostics per file

  % ---------- helpers ----------
  function H = read_gamry_header_quants(filepath, imp)
    if nargin<2 || ~isstruct(imp), imp = struct(); end
    if ~isfield(imp,'decimal_comma_to_dot'), imp.decimal_comma_to_dot = true; end
    H = struct();
    fid = fopen(filepath,'r','n');
    if fid < 0, error('Cannot open %s', filepath); end
    while true
      ln = fgetl(fid);
      if ~ischar(ln), break; end
      if imp.decimal_comma_to_dot, ln = strrep(ln, ',', '.'); end
      if ~isempty(strfind(ln, 'TABLE')), break; end
      tokens = regexp(ln, '^([A-Z0-9_]+)\s+QUANT\s+([-\+\d\.Ee]+)\s', 'tokens');
      if ~isempty(tokens)
        key = tokens{1}{1}; val = str2double(tokens{1}{2}); H.(key) = val;
      end
    end
    fclose(fid);
  end

  function [k_split, notes] = detect_split_idx(I, medwin, edge_guard_frac)
    % robust split: largest |diff| of a smoothed current, ignoring edges
    I = I(:);
    n = numel(I);
    if n < 4, k_split = max(1, floor(n/2)); notes = 'short'; return; end
    w = max(1, min(31, round(medwin)));
    Is = movmedian_oct(I, w);
    dI = diff(Is);
    guard = max(1, round(edge_guard_frac*n));
    lo = 1+guard; hi = n-1-guard;   % diff has length n-1, pick index in [lo, hi]
    if lo > hi, lo = 2; hi = n-2; end
    [~,krel] = max(abs(dI(lo:hi)));
    k = lo + krel - 1;
    % k is the index of diff, so split index is k (step1 ends at k, step2 begins at k+1)
    k_split = max(1, min(n-1, k));
    notes = '';
  end

  function y = movmedian_oct(x, w)
    % Octave-safe moving median (odd window)
    w = max(1, 2*floor((w-1)/2)+1);
    n = numel(x);
    if w == 1 || n == 0, y = x; return; end
    y = zeros(size(x));
    h = (w-1)/2;
    for i=1:n
      i1 = max(1, i-h); i2 = min(n, i+h);
      y(i) = median(x(i1:i2));
    end
  end

  function v = val_or_nan(S, field)
    if isfield(S,field) && ~isempty(S.(field)) && isfinite(S.(field)), v = S.(field);
    else v = NaN; end
  end

  function idx = pick_alpha_idx(seg_idx, alpha, mode)
    % seg_idx: vector of absolute row indices belonging to the step
    % alpha in [0,1]; alpha==0 -> FIRST sample, alpha==1 -> LAST sample (discrete)
    n = numel(seg_idx);
    if n<=0, idx = NaN; return; end
    if alpha <= 0, idx = seg_idx(1); return; end
    if alpha >= 1, idx = seg_idx(end); return; end
    pos = 1 + alpha*(n-1); % 1..n
    switch lower(mode)
      case 'floor',   ii = floor(pos);
      case 'ceil',    ii = ceil(pos);
      otherwise,      ii = round(pos);
    end
    ii = max(1, min(n, ii));
    idx = seg_idx(ii);
  end

  function process_dir(Flst, isAsc)
    if isempty(Flst), return; end
    for ff = 1:numel(Flst)
      progress(sprintf('%s: %s', Flst(ff).dir, Flst(ff).name), 0.10 + 0.70*(ff/numel(Flst)), opts);

      [t,V,I] = import_gamry_dta(Flst(ff).full, opts.import);
      H = read_gamry_header_quants(Flst(ff).full, opts.import);

      % store representative header from the first file we see
      if isempty(meta_rep)
        meta_rep = struct();
        meta_rep.SAMPLETIME = val_or_nan(H,'SAMPLETIME');
        meta_rep.TSTEP1     = val_or_nan(H,'TSTEP1');
        meta_rep.TSTEP2     = val_or_nan(H,'TSTEP2');
        meta_rep.ISTEP1     = val_or_nan(H,'ISTEP1');
        meta_rep.ISTEP2     = val_or_nan(H,'ISTEP2');
        meta_rep.DELTAI     = meta_rep.ISTEP2 - meta_rep.ISTEP1;
      end

      % local arrays
      t = t(:); V = V(:); I = I(:);
      tl = t - t(1);
      n = numel(I);

      % --- Build GLOBAL time for this direction (concatenate files) ---
      if isAsc
        tg = asc_offset + tl;  asc_offset = tg(end);
        Asc_tglob = [Asc_tglob; tg];     %#ok<AGROW>
        Asc_Vglob = [Asc_Vglob; V(:)];   %#ok<AGROW>
        Asc_Iglob = [Asc_Iglob; I(:)];   %#ok<AGROW>
      else
        tg = dsc_offset + tl;  dsc_offset = tg(end);
        Dsc_tglob = [Dsc_tglob; tg];     %#ok<AGROW>
        Dsc_Vglob = [Dsc_Vglob; V(:)];   %#ok<AGROW>
        Dsc_Iglob = [Dsc_Iglob; I(:)];   %#ok<AGROW>
      end

      % detect split index using I(t)
      [k_split, note] = detect_split_idx(I, opts.detect.medwin, opts.detect.edge_guard_frac);
      % define steps by index
      idx1 = (1:k_split).';         % step 1 samples
      idx2 = (k_split+1:n).';       % step 2 samples
      if isempty(idx1) || isempty(idx2)
        % fallback: split in the middle
        k_split = floor(n/2);
        idx1 = (1:k_split).'; idx2 = (k_split+1:n).'; note = [note '|fallback_mid'];
      end

      % build "all" rows (time-local already computed as tl)
      rows_all = [ ...
        nan(n,1), ...
        Flst(ff).fileidx * ones(n,1), ...
        [ones(numel(idx1),1); 2*ones(numel(idx2),1)], ...
        tl, V, I ...
      ];
      names_all = repmat({Flst(ff).name}, n, 1);

      % last-point rows (strictly last index of each step)
      i1_last = idx1(end); i2_last = idx2(end);
      rows_last = [ ...
        Flst(ff).fileidx, 1, tl(i1_last), V(i1_last), I(i1_last); ...
        Flst(ff).fileidx, 2, tl(i2_last), V(i2_last), I(i2_last) ...
      ];

      % alpha-picked rows (discrete at 0/1)
      i1a = pick_alpha_idx(idx1, opts.alpha, opts.alpha_round);
      i2a = pick_alpha_idx(idx2, opts.alpha, opts.alpha_round);
      rows_alpha = [ ...
        Flst(ff).fileidx, 1, tl(i1a), V(i1a), I(i1a); ...
        Flst(ff).fileidx, 2, tl(i2a), V(i2a), I(i2a) ...
      ];

      % push to direction accumulators
      if isAsc
        Asc_all      = [Asc_all; rows_all];            %#ok<AGROW>
        Asc_names    = [Asc_names; names_all];         %#ok<AGROW>
        Asc_last     = [Asc_last; rows_last];          %#ok<AGROW>
        Asc_last_names = [Asc_last_names; {Flst(ff).name; Flst(ff).name}]; %#ok<AGROW>
        Asc_alpha      = [Asc_alpha; rows_alpha];      %#ok<AGROW>
        Asc_alpha_names= [Asc_alpha_names; {Flst(ff).name; Flst(ff).name}]; %#ok<AGROW>
      else
        Dsc_all      = [Dsc_all; rows_all];            %#ok<AGROW>
        Dsc_names    = [Dsc_names; names_all];         %#ok<AGROW>
        Dsc_last     = [Dsc_last; rows_last];          %#ok<AGROW>
        Dsc_last_names = [Dsc_last_names; {Flst(ff).name; Flst(ff).name}]; %#ok<AGROW>
        Dsc_alpha      = [Dsc_alpha; rows_alpha];      %#ok<AGROW>
        Dsc_alpha_names= [Dsc_alpha_names; {Flst(ff).name; Flst(ff).name}]; %#ok<AGROW>
      end

      % QA
      Iset1 = val_or_nan(H,'ISTEP1'); Iset2 = val_or_nan(H,'ISTEP2');
      QA_rows(end+1,:) = {Flst(ff).name, 1, numel(idx1), Iset1, mean(I(idx1)), std(I(idx1)), tl(k_split), note}; %#ok<AGROW>
      QA_rows(end+1,:) = {Flst(ff).name, 2, numel(idx2), Iset2, mean(I(idx2)), std(I(idx2)), tl(k_split), note}; %#ok<AGROW>
    end
  end

  process_dir(A, true);
  process_dir(D, false);

  % fill indices
  if ~isempty(Asc_all), Asc_all(:,1) = (1:size(Asc_all,1))'; end
  if ~isempty(Dsc_all), Dsc_all(:,1) = (1:size(Dsc_all,1))'; end

  % ---------- duplicate-step check (Dsc) ----------
  if ~isempty(Dsc_last)
    Ivals = Dsc_last(:,5);                     % currents at last-point picks (two per file)
    % exclude global min and max when checking for duplicates ("besides first and last")
    if ~isempty(Ivals)
      Imin = min(Ivals); Imax = max(Ivals);
      mask_mid = (abs(Ivals - Imin) > opts.dup_tol) & (abs(Ivals - Imax) > opts.dup_tol);
      Imid = Ivals(mask_mid);
      if ~isempty(Imid)
        % detect duplicates with tolerance by rounding
        key = round(Imid/opts.dup_tol)*opts.dup_tol;
        [u,~,ic] = unique(key);
        counts = accumarray(ic, 1);
        dups = u(counts>1);
        if ~isempty(dups)
          % list unique duplicated values and how many times
          msg = sprintf('Dsc duplicate step currents detected (excluding first/last): %s', mat2str(dups));
        end
      end
    end
  end

  % ---------- build sheets ----------
  Sheets = {};

  % Meta sheet
  meta_hdr = {'SAMPLETIME_s','TSTEP1_s','TSTEP2_s','ISTEP1_A','ISTEP2_A','DeltaI_A'};
  meta_row = { meta_rep.SAMPLETIME, meta_rep.TSTEP1, meta_rep.TSTEP2, meta_rep.ISTEP1, meta_rep.ISTEP2, meta_rep.DELTAI };
  Sheets{end+1} = {'Meta', meta_hdr, meta_row};

  % QA sheet
  if ~isempty(QA_rows)
    QA_hdr = {'filename','step','N_samples','I_set_A','I_mean_A','I_std_A','t_split_s','notes'};
    Sheets{end+1} = {'Per_file_QA', QA_hdr, QA_rows};
  end

  % All samples
  headers_all = {'index','file_index','step_in_file','t_local_s','V_V','I_A','filename'};
  if ~isempty(Asc_all)
    Sheets{end+1} = {'All_Asc', headers_all, [num2cell(Asc_all) Asc_names]}; %#ok<AGROW>
  end
  if ~isempty(Dsc_all)
    Sheets{end+1} = {'All_Dsc', headers_all, [num2cell(Dsc_all) Dsc_names]}; %#ok<AGROW>
  end

  % Last-point IV (two rows per file)
  headers_pick = {'file_index','step_in_file','t_pick_s','V_V','I_A','filename'};
  if ~isempty(Asc_last)
    [~,oa] = sortrows(Asc_last,[1 2]);
    Sheets{end+1} = {'IV_last_Asc', headers_pick, [num2cell(Asc_last(oa,:)) Asc_last_names(oa)]}; %#ok<AGROW>
  end
  if ~isempty(Dsc_last)
    [~,od] = sortrows(Dsc_last,[1 2]);
    Sheets{end+1} = {'IV_last_Dsc', headers_pick, [num2cell(Dsc_last(od,:)) Dsc_last_names(od)]}; %#ok<AGROW>
  end

  % Alpha-point IV (two rows per file)
  if ~isempty(Asc_alpha)
    [~,oa] = sortrows(Asc_alpha,[1 2]);
    Sheets{end+1} = {'IV_alpha_Asc', headers_pick, [num2cell(Asc_alpha(oa,:)) Asc_alpha_names(oa)]}; %#ok<AGROW>
  end
  if ~isempty(Dsc_alpha)
    [~,od] = sortrows(Dsc_alpha,[1 2]);
    Sheets{end+1} = {'IV_alpha_Dsc', headers_pick, [num2cell(Dsc_alpha(od,:)) Dsc_alpha_names(od)]}; %#ok<AGROW>
  end

  % ---------- plots ----------
  progress('Plotting...', 0.88, opts);

  % Full dataset: V vs I (Asc & Dsc)
  hVI = figure('Name','Polarización: V vs I (full samples)'); hold on; grid on;
  legend_str = {};
  if ~isempty(Asc_all)
    plot(Asc_all(:,6), Asc_all(:,5), '.', 'MarkerSize', 6);
    legend_str{end+1}='Asc (all)';
  end
  if ~isempty(Dsc_all)
    plot(Dsc_all(:,6), Dsc_all(:,5), '.', 'MarkerSize', 6);
    legend_str{end+1}='Dsc (all)';
  end
  xlabel('Current (A)'); ylabel('Voltage (V)');
  if ~isempty(legend_str), legend(legend_str); end
  title('Polarización: V vs I (muestras completas)');

  % Last-point IV: exactly two points per file & direction
  hIVlast = figure('Name','Polarización: IV (último punto por escalón)'); hold on; grid on;
  legend_str = {};
  if ~isempty(Asc_last)
    [Iasc, ia] = sort(Asc_last(:,5)); Vasc = Asc_last(ia,4);
    plot(Iasc, Vasc, '-o', 'LineWidth', 1.0, 'MarkerSize', 3);
    legend_str{end+1}='Asc (last points)';
  end
  if ~isempty(Dsc_last)
    [Idsc, id] = sort(Dsc_last(:,5)); Vdsc = Dsc_last(id,4);
    plot(Idsc, Vdsc, '-o', 'LineWidth', 1.0, 'MarkerSize', 3);
    legend_str{end+1}='Dsc (last points)';
  end
  xlabel('Current (A)'); ylabel('Voltage (V)');
  if ~isempty(legend_str), legend(legend_str); end
  title('Polarización: IV (dos puntos por archivo: fin de escalón 1 y 2)');

  % ---- Dual-axis global plots per direction (V & I vs t) ----
  if ~isempty(Asc_tglob)
    hAscDual = figure('Name','Polarización Asc: V and I vs t (global)');
    [haxA, hVA, hIA] = plotyy(Asc_tglob, Asc_Vglob, Asc_tglob, Asc_Iglob);
    grid(haxA(1),'on'); grid(haxA(2),'on');
    xlabel('Time (s)');
    ylabel(haxA(1), 'Voltage (V)');
    ylabel(haxA(2), 'Current (A)');
    title('Polarización Asc: V and I vs t (global)');
    legend([hVA; hIA], {'V (left)','I (right)'});
  end
  if ~isempty(Dsc_tglob)
    hDscDual = figure('Name','Polarización Dsc: V and I vs t (global)');
    [haxD, hVD, hID] = plotyy(Dsc_tglob, Dsc_Vglob, Dsc_tglob, Dsc_Iglob);
    grid(haxD(1),'on'); grid(haxD(2),'on');
    xlabel('Time (s)');
    ylabel(haxD(1), 'Voltage (V)');
    ylabel(haxD(2), 'Current (A)');
    title('Polarización Dsc: V and I vs t (global)');
    legend([hVD; hID], {'V (left)','I (right)'});
  end

  % ---------- Impedancia diferencial Zdiff = dV/dI (desde puntos finales) ----------
  % Usa Asc_last y Dsc_last. Cada fila: [fileidx, step_in_file, t_last_s, V_last, I_last]
  hZdiff = figure('Name','Polarización: Z_{diff} = dV/dI'); hold on; grid on;
  legZd = {};
  % Ventana de suavizado (impar)
  w = opts.diff_smooth_pts; if isempty(w), w = 1; end
  if mod(w,2)==0, w = w+1; end
  % Asc
  if ~isempty(Asc_last)
    [Iasc, ia] = sort(Asc_last(:,5)); Vasc = Asc_last(ia,4);
    Vs = Vasc; if w>1, Vs = smooth_ma(Vasc, w); end
    if numel(Iasc) >= 2
      dI = diff(Iasc); dV = diff(Vs);
      m = abs(dI) > eps;
      I_mid = (Iasc(1:end-1)+Iasc(2:end))/2;
      plot(I_mid(m), dV(m)./dI(m), '-o', 'LineWidth', 1.0, 'MarkerSize', 4);
      legZd{end+1} = 'Asc dV/dI';
    end
  end
  % Dsc
  if ~isempty(Dsc_last)
    [Idsc, id] = sort(Dsc_last(:,5)); Vdsc = Dsc_last(id,4);
    Vs = Vdsc; if w>1, Vs = smooth_ma(Vdsc, w); end
    if numel(Idsc) >= 2
      dI = diff(Idsc); dV = diff(Vs);
      m = abs(dI) > eps;
      I_mid = (Idsc(1:end-1)+Idsc(2:end))/2;
      plot(I_mid(m), dV(m)./dI(m), '-o', 'LineWidth', 1.0, 'MarkerSize', 4);
      legZd{end+1} = 'Dsc dV/dI';
    end
  end
  xlabel('Current (A)'); ylabel('dV/dI (\Omega)');
  title('Polarización: impedancia diferencial');
  if ~isempty(legZd), legend(legZd); end

  % ---------- export ----------
  R = struct('default_filename', opts.export_filename, 'sheets', {Sheets});
  if isfield(opts,'do_export') && opts.do_export
    progress('Exporting workbook...', 0.94, opts);
    export_workbook(R, opts.export_filename, opts.progress_cb);
  end

  % ---------- SVG export ----------
  if isfield(opts,'save_svg') && opts.save_svg
    try
      if ~exist(opts.svg_dir,'dir'), mkdir(opts.svg_dir); end
      % Set painters to reduce gl2ps issues and use safe exporter with fallback
      if exist('hVI','var'),      export_svg_safe(hVI,     fullfile(opts.svg_dir, sprintf('POL_V_vs_I_full_%s.svg', stamp)), 'Width', 1600, 'Height', 1000); end
      if exist('hIVlast','var'),  export_svg_safe(hIVlast, fullfile(opts.svg_dir, sprintf('POL_IV_last_%s.svg', stamp)),     'Width', 1600, 'Height', 1000); end
      if exist('hAscDual','var'), export_svg_safe(hAscDual, fullfile(opts.svg_dir, sprintf('POL_Asc_VI_vs_t_%s.svg', stamp)), 'Width', 1600, 'Height', 1000); end
      if exist('hDscDual','var'), export_svg_safe(hDscDual, fullfile(opts.svg_dir, sprintf('POL_Dsc_VI_vs_t_%s.svg', stamp)), 'Width', 1600, 'Height', 1000); end
      if exist('hZdiff','var'),   export_svg_safe(hZdiff,  fullfile(opts.svg_dir, sprintf('POL_Zdiff_%s.svg', stamp)),        'Width', 1600, 'Height', 1000); end
    catch ME
      warning('SVG export failed: %s', ME.message);
    end
  end

  progress('Done', 1.00, opts);
end

function export_svg_safe(fig, outfn, varargin)
  % Robust SVG export: try painters; on failure, temporarily switch to gnuplot.
  try
    if ishghandle(fig)
      try, set(fig, 'renderer', 'painters'); catch, end
      drawnow();
    end
    export_svg_scaled(fig, outfn, varargin{:});
  catch
    try
      oldtk = graphics_toolkit();
      graphics_toolkit('gnuplot');
      drawnow();
      export_svg_scaled(fig, outfn, varargin{:});
      graphics_toolkit(oldtk);
    catch ME
      warning('SVG export failed (safe): %s', ME.message);
    end
  end
end

function y = smooth_ma(y, w)
  % Moving-average smoothing with odd window length w >= 1.
  if nargin<2 || isempty(w), w = 1; end
  w = max(1, round(w));
  if mod(w,2)==0, w = w+1; end
  if w<=1 || numel(y)<3, return; end
  k = (w-1)/2;
  ypad = [repmat(y(1),k,1); y(:); repmat(y(end),k,1)];
  ker = ones(w,1)/w;
  y = conv(ypad, ker, 'valid');
end


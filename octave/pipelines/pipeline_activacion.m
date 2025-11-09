function R = pipeline_activacion(folder, opts)
% PIPELINE_ACTIVACION
% Processes activation runs:
%   Activacion_(Asc|Dsc)_(proto)_#CYCLE_#FILE.DTA
% Each file contains TWO current steps. We:
%   - Stitch Asc then Dsc within each CYCLE into a single global timeline
%   - Compute per-sample step_in_file (1 or 2) using detect_plateaus()
%   - Plot global V vs t, global I vs t
%   - Plot cycle overlays (V vs t) for Asc and Dsc (local time in each cycle)
%   - Export (last step): All + one sheet per (Cycle, Direction)
%   - Alpha point per step (default last point) and export summaries,
%     plus plots (V vs t) with small dots joined by lines.
%   - NEW: Extract header meta from FIRST file (IPRESTEP, TPRESTEP, ISTEP1,
%     TSTEP1, ISTEP2, TSTEP2, SAMPLETIME) into a 'Meta' sheet.
%   - NEW: Add 'SummaryCycles' sheet with number of cycles processed.
%   - NEW: Export SVGs for all figures via export_svg_scaled().

  addpath(fullfile(pwd,'lib'));
  clear i j  % Octave safety

  % ---------- defaults ----------
  stamp = datestr(now, 'yyyymmdd_HHMM');
  defaults = struct( ...
    'alpha', 1.0, ...
    'do_export', true, ...
    'export_filename', sprintf('Gamry_Activacion_%s.xlsx', stamp), ...
    'dir_priority', {{'Asc','Dsc'}}, ...
    'progress_cb', [], ...
    'save_svg', true, ...
    'svg_dir', fullfile(folder, 'plots_svg') ...
  );
  defaults.import = struct('header_guess_lines',160,'min_numeric_cols',5,'decimal_comma_to_dot',true);
  defaults.detect = struct('smooth_frac',0.01,'edge_guard',0.05,'min_len_frac',0.05);
  if nargin < 2 || ~isstruct(opts) || isempty(opts), opts = struct(); end
  opts = merge_opts(defaults, opts);

  % clamp alpha [0,1]
  if ~isfield(opts,'alpha') || isempty(opts.alpha), opts.alpha = 1.0; end
  opts.alpha = max(0,min(1,opts.alpha));

  % ---------- discover files ----------
  progress('Scanning files...', 0.04, opts);
  pat = '^Activacion_(Asc|Dsc)_(.+?)_#(\d+)_#(\d+)\.DTA$';
  F = list_files(folder, pat);
  if isempty(F), error('No activation files found in: %s', folder); end

  % parse table
  T = struct('name',{},'full',{},'dir',{},'proto',{},'cycle',{},'fileidx',{});
  for kk = 1:numel(F)
    tok = F(kk).tokens;  % {dir, proto, cycle, fileidx}
    T(end+1).name    = F(kk).name;      %#ok<AGROW>
    T(end).full      = F(kk).full;
    T(end).dir       = tok{1};
    T(end).proto     = tok{2};
    T(end).cycle     = str2double(tok{3});
    T(end).fileidx   = str2double(tok{4});
  end

  cycles = sort(unique([T.cycle]));

  % ---------- FIRST-FILE META (from the first file found) ----------
  MetaSheet = {};
  try
    H = read_activation_header_meta(T(1).full, opts.import);
    meta_hdr = {'IPRESTEP_A','TPRESTEP_s','ISTEP1_A','TSTEP1_s','ISTEP2_A','TSTEP2_s','SAMPLETIME_s', ...
                'FRAMEWORKVERSION','INSTRUMENTVERSION','PSTATSERIALNO','DATE','TIME','PSTAT'};
    meta_row = {H.IPRESTEP, H.TPRESTEP, H.ISTEP1, H.TSTEP1, H.ISTEP2, H.TSTEP2, H.SAMPLETIME, ...
                H.FRAMEWORKVERSION, H.INSTRUMENTVERSION, H.PSTATSERIALNO, H.DATE, H.TIME, H.PSTAT};
    MetaSheet = { {'Meta', meta_hdr, meta_row} };
  catch ME
    warning('Meta extraction failed: %s', ME.message);
  end

  % ---------- accumulation (global + per-cycle sheets) ----------
  progress('Reading & stitching...', 0.10, opts);

  % global accumulators
  G_t_global = []; G_t_local = []; G_V = []; G_I = []; G_step = []; G_fileidx = []; G_cycle = [];
  global_offset = 0;

  % per-cycle sheets
  Sheets = {};
  headers = {'index','file_index','cycle_index','step_in_file','t_global_s','t_local_s','V_V','I_A'};

  % overlays
  AscOverlays = struct('tcells',{{}},'vcells',{{}},'labels',{{}});
  DscOverlays = struct('tcells',{{}},'vcells',{{}},'labels',{{}});

  % alpha-summary accumulators
  % columns: cycle, file_index, step_in_file(1/2), t_cycle_s, V_V, I_A, filename
  SumAsc = []; SumAsc_names = {};
  SumDsc = []; SumDsc_names = {};

  % loop cycles
  for cc = 1:numel(cycles)
    cyc = cycles(cc);

    A = T(strcmp({T.dir},'Asc') & [T.cycle]==cyc);
    D = T(strcmp({T.dir},'Dsc') & [T.cycle]==cyc);
    [~,aord] = sort([A.fileidx]); A = A(aord);
    [~,dord] = sort([D.fileidx]); D = D(dord);

    cyc_list = {};
    for dd = 1:numel(opts.dir_priority)
      dtag = opts.dir_priority{dd};
      if strcmp(dtag,'Asc'), cyc_list{end+1} = A; %#ok<AGROW>
      elseif strcmp(dtag,'Dsc'), cyc_list{end+1} = D; %#ok<AGROW>
      end
    end

    CycAsc_rows = [];  % [idx, fileidx, cycle, step, tglob, tloc, V, I]
    CycDsc_rows = [];

    asc_local_t = []; asc_local_v = [];
    dsc_local_t = []; dsc_local_v = [];
    asc_accum_offset = 0; dsc_accum_offset = 0;

    for dl = 1:numel(cyc_list)
      Fdir = cyc_list{dl};
      if isempty(Fdir)
        progress(sprintf('Cycle %d: %s missing (skipping)',cyc, opts.dir_priority{dl}), 0.10, opts);
        continue;
      end

      for ff = 1:numel(Fdir)
        progress(sprintf('Cycle %d %s: %s',cyc, Fdir(ff).dir, Fdir(ff).name), ...
                 0.12 + 0.60*((cc-1 + ff/max(1,numel(Fdir)))/max(1,numel(cycles))), opts);

        [t,V,I] = import_gamry_dta(Fdir(ff).full, opts.import);
        tl = t(:) - t(1);

        % detect 2 plateaus in this file to tag each sample with step 1/2
        segs = detect_plateaus(I(:), 2, opts.detect);
        % fallback boundary if detection under-delivers
        if size(segs,1) < 2
          hb = max(2, round(numel(tl)/2));
          segs = [1 hb; hb+1 numel(tl)];
        end

        % build sample step vector (1/2)
        stepvec = ones(size(tl));
        step_boundary = segs(1,2);
        stepvec((step_boundary+1):end) = 2;

        % append to GLOBAL arrays (continuous time)
        tg = global_offset + tl;
        G_t_global = [G_t_global; tg];                  %#ok<AGROW>
        G_t_local  = [G_t_local;  tl];                  %#ok<AGROW>
        G_V        = [G_V; V(:)];                       %#ok<AGROW>
        G_I        = [G_I; I(:)];                       %#ok<AGROW>
        G_step     = [G_step; stepvec(:)];              %#ok<AGROW>
        G_fileidx  = [G_fileidx; Fdir(ff).fileidx*ones(numel(tl),1)]; %#ok<AGROW>
        G_cycle    = [G_cycle; cyc*ones(numel(tl),1)];  %#ok<AGROW>
        global_offset = tg(end);

        % rows for per-cycle sheets
        rows = [ ...
          nan(numel(tl),1), ...
          Fdir(ff).fileidx * ones(numel(tl),1), ...
          cyc * ones(numel(tl),1), ...
          stepvec(:), ...
          tg(:), ...
          tl(:), ...
          V(:), ...
          I(:) ...
        ];

        % cycle-local concatenation for overlays and for alpha timing
        if strcmp(Fdir(ff).dir,'Asc')
          tcl = asc_accum_offset + tl;
          asc_local_t = [asc_local_t; tcl]; %#ok<AGROW>
          asc_local_v = [asc_local_v; V(:)]; %#ok<AGROW>
          asc_accum_offset = tcl(end);
          CycAsc_rows = [CycAsc_rows; rows]; %#ok<AGROW>
        else
          tcl = dsc_accum_offset + tl;
          dsc_local_t = [dsc_local_t; tcl]; %#ok<AGROW>
          dsc_local_v = [dsc_local_v; V(:)]; %#ok<AGROW>
          dsc_accum_offset = tcl(end);
          CycDsc_rows = [CycDsc_rows; rows]; %#ok<AGROW>
        end

        % alpha points per step (cycle-local time)
        alpha_idx = @(seg) max(seg(1), min(seg(2), seg(1) + round(opts.alpha * max(0, (seg(2)-seg(1))))));
        i1 = alpha_idx(segs(1,:));
        i2 = alpha_idx(segs(2,:));

        if strcmp(Fdir(ff).dir,'Asc')
          SumAsc(end+1,:) = [cyc, Fdir(ff).fileidx, 1, tcl(i1), V(i1), I(i1)]; %#ok<AGROW>
          SumAsc_names{end+1,1} = Fdir(ff).name; %#ok<AGROW>
          SumAsc(end+1,:) = [cyc, Fdir(ff).fileidx, 2, tcl(i2), V(i2), I(i2)]; %#ok<AGROW>
          SumAsc_names{end+1,1} = Fdir(ff).name; %#ok<AGROW>
        else
          SumDsc(end+1,:) = [cyc, Fdir(ff).fileidx, 1, tcl(i1), V(i1), I(i1)]; %#ok<AGROW>
          SumDsc_names{end+1,1} = Fdir(ff).name; %#ok<AGROW>
          SumDsc(end+1,:) = [cyc, Fdir(ff).fileidx, 2, tcl(i2), V(i2), I(i2)]; %#ok<AGROW>
          SumDsc_names{end+1,1} = Fdir(ff).name; %#ok<AGROW>
        end
      end % files
    end % direction within cycle

    % finalize per-cycle sheets
    if ~isempty(CycAsc_rows)
      CycAsc_rows(:,1) = (1:size(CycAsc_rows,1))';
      Sheets{end+1} = {sprintf('Cycle%02d_Asc',cyc), headers, num2cell(CycAsc_rows)}; %#ok<AGROW>
      AscOverlays.tcells{end+1} = asc_local_t; %#ok<AGROW>
      AscOverlays.vcells{end+1} = asc_local_v; %#ok<AGROW>
      AscOverlays.labels{end+1} = sprintf('Cycle %d', cyc); %#ok<AGROW>
    end
    if ~isempty(CycDsc_rows)
      CycDsc_rows(:,1) = (1:size(CycDsc_rows,1))';
      Sheets{end+1} = {sprintf('Cycle%02d_Dsc',cyc), headers, num2cell(CycDsc_rows)}; %#ok<AGROW>
      DscOverlays.tcells{end+1} = dsc_local_t; %#ok<AGROW>
      DscOverlays.vcells{end+1} = dsc_local_v; %#ok<AGROW>
      DscOverlays.labels{end+1} = sprintf('Cycle %d', cyc); %#ok<AGROW>
    end
  end % cycles

  % ---------- build "All" sheet (global joined order) ----------
  N = numel(G_t_global);
  All_rows = [ (1:N)', G_fileidx, G_cycle, G_step, G_t_global, G_t_local, G_V, G_I ];
  Sheets = [ { {'All', headers, num2cell(All_rows)} }, Sheets ];

  % ---------- Summary sheets ----------
  % Alpha summaries
  if ~isempty(SumAsc)
    [~,ordA] = sortrows(SumAsc, [1 2 3]);
    SumAsc = SumAsc(ordA,:); SumAsc_names = SumAsc_names(ordA);
    headersS = {'cycle_index','file_index','step_in_file','t_cycle_s','V_V','I_A','filename'};
    Sheets{end+1} = {'Summary_Asc', headersS, [num2cell(SumAsc) SumAsc_names]}; %#ok<AGROW>
  end
  if ~isempty(SumDsc)
    [~,ordD] = sortrows(SumDsc, [1 2 3]);
    SumDsc = SumDsc(ordD,:); SumDsc_names = SumDsc_names(ordD);
    headersS = {'cycle_index','file_index','step_in_file','t_cycle_s','V_V','I_A','filename'};
    Sheets{end+1} = {'Summary_Dsc', headersS, [num2cell(SumDsc) SumDsc_names]}; %#ok<AGROW>
  end
  % Count of cycles processed
  try
    ncycles = numel(cycles);
    Sheets{end+1} = {'SummaryCycles', {'n_cycles','first_cycle','last_cycle'}, {ncycles, cycles(1), cycles(end)}};
  catch
    Sheets{end+1} = {'SummaryCycles', {'n_cycles'}, {numel(cycles)}};
  end

  % ---------- prepend Meta if available ----------
  if ~isempty(MetaSheet)
    Sheets = [ MetaSheet, Sheets ];
  end

  % ---------- plots ----------
  progress('Plotting...', 0.86, opts);

  % Global V vs t
  hVt = figure('Name','Activacion: V vs t (global joined)');
  plot(G_t_global, G_V, 'LineWidth', 1.2);
  xlabel('Time (s)'); ylabel('Voltage (V)'); grid on;
  title('Activación: V vs t (Asc→Dsc joined per cycle)');

  % Global I vs t
  hIt = figure('Name','Activacion: I vs t (global joined)');
  plot(G_t_global, G_I, 'LineWidth', 1.2);
  xlabel('Time (s)'); ylabel('Current (A)'); grid on;
  title('Activación: I vs t (Asc→Dsc joined per cycle)');

  % Overlays: Asc (V vs t_cycle)
  hAsc = [];
  if ~isempty(AscOverlays.tcells)
    hAsc = figure('Name','Activacion: Overlays Asc (V vs t_{cycle})'); hold on;
    for kk = 1:numel(AscOverlays.tcells)
      plot(AscOverlays.tcells{kk}, AscOverlays.vcells{kk}, 'LineWidth', 1.2);
    end
    xlabel('t_{cycle} (s)'); ylabel('Voltage (V)'); grid on; hold off;
    legend(AscOverlays.labels);
    title('Activación: Overlays Asc (V vs t_{cycle})');
  end

  % Overlays: Dsc (V vs t_cycle)
  hDsc = [];
  if ~isempty(DscOverlays.tcells)
    hDsc = figure('Name','Activacion: Overlays Dsc (V vs t_{cycle})'); hold on;
    for kk = 1:numel(DscOverlays.tcells)
      plot(DscOverlays.tcells{kk}, DscOverlays.vcells{kk}, 'LineWidth', 1.2);
    end
    xlabel('t_{cycle} (s)'); ylabel('Voltage (V)'); grid on; hold off;
    legend(DscOverlays.labels);
    title('Activación: Overlays Dsc (V vs t_{cycle})');
  end

  % Summary plots (alpha points)
  hAscS = [];
  if ~isempty(SumAsc)
    hAscS = figure('Name','Activacion: Summary Asc (alpha points: V vs t_{cycle})'); hold on; grid on;
    cyclesA = unique(SumAsc(:,1));
    for c = cyclesA(:)'
      m = SumAsc(:,1) == c;
      [~,ord] = sortrows(SumAsc(m,[1 2 3]));
      tc = SumAsc(m,4); Vc = SumAsc(m,5);
      tc = tc(ord);     Vc = Vc(ord);
      plot(tc, Vc, '-o', 'LineWidth', 1.0, 'MarkerSize', 4);
    end
    xlabel('t_{cycle} (s)'); ylabel('Voltage (V)');
    title(sprintf('Activación Asc summary (alpha=%.2f)', opts.alpha));
    hold off;
  end

  hDscS = [];
  if ~isempty(SumDsc)
    hDscS = figure('Name','Activacion: Summary Dsc (alpha points: V vs t_{cycle})'); hold on; grid on;
    cyclesD = unique(SumDsc(:,1));
    for c = cyclesD(:)'
      m = SumDsc(:,1) == c;
      [~,ord] = sortrows(SumDsc(m,[1 2 3]));
      tc = SumDsc(m,4); Vc = SumDsc(m,5);
      tc = tc(ord);     Vc = Vc(ord);
      plot(tc, Vc, '-o', 'LineWidth', 1.0, 'MarkerSize', 4);
    end
    xlabel('t_{cycle} (s)'); ylabel('Voltage (V)');
    title(sprintf('Activación Dsc summary (alpha=%.2f)', opts.alpha));
    hold off;
  end

  % ---------- export LAST (optional) ----------
  R = struct('default_filename', opts.export_filename, 'sheets', {Sheets});
  if isfield(opts,'do_export') && opts.do_export
    progress('Exporting workbook...', 0.92, opts);
    export_workbook(R, opts.export_filename, opts.progress_cb);
  end

  % ---------- SVG export ----------
  if opts.save_svg
    try
      if ~exist(opts.svg_dir,'dir'), mkdir(opts.svg_dir); end
      export_svg_scaled(hVt,  fullfile(opts.svg_dir, sprintf('Activacion_V_vs_t_%s.svg', stamp)));
      export_svg_scaled(hIt,  fullfile(opts.svg_dir, sprintf('Activacion_I_vs_t_%s.svg', stamp)));
      if ~isempty(hAsc) && ishandle(hAsc)
        export_svg_scaled(hAsc, fullfile(opts.svg_dir, sprintf('Activacion_Overlays_Asc_%s.svg', stamp)));
      end
      if ~isempty(hDsc) && ishandle(hDsc)
        export_svg_scaled(hDsc, fullfile(opts.svg_dir, sprintf('Activacion_Overlays_Dsc_%s.svg', stamp)));
      end
      if ~isempty(hAscS) && ishandle(hAscS)
        export_svg_scaled(hAscS, fullfile(opts.svg_dir, sprintf('Activacion_Summary_Asc_%s.svg', stamp)));
      end
      if ~isempty(hDscS) && ishandle(hDscS)
        export_svg_scaled(hDscS, fullfile(opts.svg_dir, sprintf('Activacion_Summary_Dsc_%s.svg', stamp)));
      end
    catch ME
      warning('SVG export failed: %s', ME.message);
    end
  end

  progress('Done', 1.00, opts);
end

% ================= helpers =================

function H = read_activation_header_meta(filepath, imp)
% Parse selected activation header fields up to CURVE table.
% Normalizes decimal commas to dots. Missing -> NaN or ''.
  if nargin < 2 || ~isstruct(imp), imp = struct(); end
  if ~isfield(imp,'decimal_comma_to_dot'), imp.decimal_comma_to_dot = true; end

  H = struct('IPRESTEP',NaN,'TPRESTEP',NaN,'ISTEP1',NaN,'TSTEP1',NaN, ...
             'ISTEP2',NaN,'TSTEP2',NaN,'SAMPLETIME',NaN, ...
             'FRAMEWORKVERSION','', 'INSTRUMENTVERSION','', 'PSTATSERIALNO','', ...
             'DATE','', 'TIME','', 'PSTAT','');

  fid = fopen(filepath,'r','n');
  if fid < 0, error('Cannot open %s', filepath); end
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if imp.decimal_comma_to_dot, ln = strrep(ln, ',', '.'); end
    if ~isempty(strfind(ln, 'CURVE')), break; end % stop before table

    % QUANTs
    tokQ = regexp(ln, '^([A-Z0-9_]+)\s+QUANT\s+([-+\d\.Ee]+)\s', 'tokens');
    if ~isempty(tokQ)
      key = tokQ{1}{1}; val = str2double(tokQ{1}{2});
      switch key
        case 'IPRESTEP',    H.IPRESTEP   = val;
        case 'TPRESTEP',    H.TPRESTEP   = val;
        case 'ISTEP1',      H.ISTEP1     = val;
        case 'TSTEP1',      H.TSTEP1     = val;
        case 'ISTEP2',      H.ISTEP2     = val;
        case 'TSTEP2',      H.TSTEP2     = val;
        case 'SAMPLETIME',  H.SAMPLETIME = val;
      end
      continue;
    end

    % LABELs of interest (versions/serial/date/time/title/pstat)
    tokL = regexp(ln, '^([A-Z0-9_]+)\s+LABEL\s+(.+?)\s*$', 'tokens');
    if ~isempty(tokL)
      key = tokL{1}{1}; sval = strtrim(tokL{1}{2});
      switch key
        case 'FRAMEWORKVERSION',   H.FRAMEWORKVERSION  = sval;
        case 'INSTRUMENTVERSION',  H.INSTRUMENTVERSION = sval;
        case 'PSTATSERIALNO',      H.PSTATSERIALNO     = sval;
        case 'DATE',               H.DATE              = sval;
        case 'TIME',               H.TIME              = sval;
        case 'PSTAT',              H.PSTAT             = sval;
      end
      continue;
    end
  end
  fclose(fid);
end


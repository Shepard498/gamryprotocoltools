function pipeline_activation(dirs, opts)
  if nargin < 2, opts = struct(); end
  C = config();
  opts = merge_opts(struct('svg', true, 'png', false, 'fig_visible', C.fig_visible, 'size_px', C.fig_size_px), opts);

  files = [ dir(fullfile(dirs.raw, 'Activaci*_*.*DTA')) ];
  if isempty(files), return; end

  S = sort_tokens(files, 'activation');
  xlsx = fullfile(dirs.processed, 'Activation_results.xlsx');
  sheets = {};

  G_raw = struct('global_idx', [], 'cycle', [], 'step', [], 't_s', [], 'V_V', [], 'I_A', []);
  G_last = struct('global_idx', [], 'cycle', [], 'step', [], 't_s', [], 'V_V', [], 'I_A', []);
  asc_cycles = {}; dsc_cycles = {};

  meta_first = struct();
  area_used = C.active_area_cm2_default; area_warn = false;

  fig_global = figure('Name','ACT Global timeline'); set(fig_global,'Visible',opts.fig_visible);
  tiledlayout(fig_global,2,1);
  ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); ylabel(ax1,'V (V)');
  ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); ylabel(ax2,'I (A)'); xlabel(ax2,'t (s)');

  fig_global_last = figure('Name','ACT Global timeline (last points)'); set(fig_global_last,'Visible',opts.fig_visible);
  hold(fig_global_last,'on'); grid on; xlabel('t (s)'); ylabel('V (V)');

  asc_local_fig = figure('Name','ACT ASC overlay (local t)'); set(asc_local_fig,'Visible',opts.fig_visible); hold on; grid on; xlabel('t (s)'); ylabel('V (V)');
  dsc_local_fig = figure('Name','ACT DSC overlay (local t)'); set(dsc_local_fig,'Visible',opts.fig_visible); hold on; grid on; xlabel('t (s)'); ylabel('V (V)');

  asc_local_last = figure('Name','ACT ASC overlay (last points)'); set(asc_local_last,'Visible',opts.fig_visible); hold on; grid on; xlabel('I (A)'); ylabel('V (V)');
  dsc_local_last = figure('Name','ACT DSC overlay (last points)'); set(dsc_local_last,'Visible',opts.fig_visible); hold on; grid on; xlabel('I (A)'); ylabel('V (V)');

  gidx = 0; step_counter = 0; cur_cycle = NaN;
  for k = 1:numel(S)
    fp = S(k).path; nm = S(k).name;
    Ctab = import_gamry_curves(fp);
    M = read_dta_metadata(fp);
    if k==1, meta_first = M; end
    if isfield(M, 'area'), area_used = M.area; else, area_warn = true; end

    plot(ax1, Ctab.t_s, Ctab.V_V, '-','DisplayName', nm);
    plot(ax2, Ctab.t_s, Ctab.I_A, '-','DisplayName', nm);

    if isnan(cur_cycle), cur_cycle = S(k).cycle; end
    if S(k).cycle ~= cur_cycle
      cur_cycle = S(k).cycle;
    end
    n = numel(Ctab.t_s);
    G_raw.global_idx = [G_raw.global_idx; (gidx+1:gidx+n)'];
    step_ids = zeros(n,1);
    cyc_ids = repmat(S(k).cycle, n, 1);
    G_raw.cycle = [G_raw.cycle; cyc_ids];
    G_raw.step  = [G_raw.step; step_ids];
    G_raw.t_s   = [G_raw.t_s; Ctab.t_s];
    G_raw.V_V   = [G_raw.V_V; Ctab.V_V];
    G_raw.I_A   = [G_raw.I_A; Ctab.I_A];
    gidx = gidx + n;

    cfg = C.step;
    P = last_points(Ctab.t_s, Ctab.I_A, Ctab.V_V, cfg);
    step_counter = step_counter + numel(P.I_last_A);

    G_last.global_idx = [G_last.global_idx; (numel(G_last.global_idx)+(1:numel(P.t_last_s)))'];
    G_last.cycle = [G_last.cycle; repmat(S(k).cycle, numel(P.t_last_s),1)];
    G_last.step  = [G_last.step; P.step_idx];
    G_last.t_s   = [G_last.t_s; P.t_last_s];
    G_last.V_V   = [G_last.V_V; P.V_last_V];
    G_last.I_A   = [G_last.I_A; P.I_last_A];

    plot(fig_global_last, P.t_last_s, P.V_last_V, 'o-','DisplayName', nm);

    if S(k).dirflag
      plot(asc_local_fig, Ctab.t_s - Ctab.t_s(1), Ctab.V_V, '-','DisplayName', nm);
      plot(asc_local_last, P.I_last_A, P.V_last_V, 'o-','DisplayName', nm);
      asc_cycles{end+1} = struct('name', nm, 'P', P);
    else
      plot(dsc_local_fig, Ctab.t_s - Ctab.t_s(1), Ctab.V_V, '-','DisplayName', nm);
      plot(dsc_local_last, P.I_last_A, P.V_last_V, 'o-','DisplayName', nm);
      dsc_cycles{end+1} = struct('name', nm, 'P', P);
    end
  end

  if area_warn
    warning('AREA not found in metadata for some files. Using default %.3g cm^2', area_used);
  end

  meta_cells = {
    'DATE', safe_get(meta_first,'date',''); 'TIME', safe_get(meta_first,'time','');
    'ISTEP1', safe_get(meta_first,'istep1', NaN); 'TSTEP1', safe_get(meta_first,'tstep1', NaN);
    'ISTEP2', safe_get(meta_first,'istep2', NaN); 'TSTEP2', safe_get(meta_first,'tstep2', NaN);
    'SAMPLETIME', safe_get(meta_first,'sampletime', NaN);
    'AREA_cm2', area_used;
  };
  sheets(end+1,:) = {'META', meta_cells};

  sheets(end+1,:) = {'GLOBAL_raw', table_cells(G_raw)};
  sheets(end+1,:) = {'GLOBAL_last', table_cells(G_last)};

  for c = 1:numel(asc_cycles)
    nm = sprintf('ASC_cyc%02d', c);
    P = asc_cycles{c}.P; T = struct('I_A',P.I_last_A,'V_V',P.V_last_V,'t_s',P.t_last_s);
    sheets(end+1,:) = {nm, table_cells(T)};
  end
  for c = 1:numel(dsc_cycles)
    nm = sprintf('DSC_cyc%02d', c);
    P = dsc_cycles{c}.P; T = struct('I_A',P.I_last_A,'V_V',P.V_last_V,'t_s',P.t_last_s);
    sheets(end+1,:) = {nm, table_cells(T)};
  end

  Smet = activation_metrics(G_last, area_used);
  summary_cells = {
    'n_cycles', Smet.n_cycles; 'n_steps', Smet.n_steps; 't_total_s', Smet.t_total_s;
    'I_min_A', Smet.I_range(1); 'I_max_A', Smet.I_range(2);
    'V_min_V', Smet.V_range(1); 'V_max_V', Smet.V_range(2);
    'J_min_Acm2', Smet.J_range(1); 'J_max_Acm2', Smet.J_range(2);
  };
  sheets(end+1,:) = {'SUMMARY', summary_cells};

  export_excel(xlsx, sheets);

  if opts.svg
    save_svg(fig_global, fullfile(dirs.plots,'ACT_global.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_global_last, fullfile(dirs.plots,'ACT_global_last.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(asc_local_fig, fullfile(dirs.plots,'ACT_ASC_overlay.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(dsc_local_fig, fullfile(dirs.plots,'ACT_DSC_overlay.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(asc_local_last, fullfile(dirs.plots,'ACT_ASC_last.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(dsc_local_last, fullfile(dirs.plots,'ACT_DSC_last.svg'), struct('png',opts.png,'size_px',opts.size_px));
  end
end

function v = safe_get(S, f, dv)
  if isfield(S, f), v = S.(f); else, v = dv; end
end

function cells = table_cells(S)
  f = fieldnames(S);
  hdr = f(:)';
  n = numel(S.(f{1}));
  cells = [hdr; cell(n, numel(f))];
  for j = 1:numel(f)
    v = S.(f{j});
    if iscell(v), v = v(:); else, v = num2cell(v(:)); end
    cells(2:end,j) = v;
  end
end


function pipeline_polarization(dirs, opts)
  if nargin < 2, opts = struct(); end
  C = config();
  opts = merge_opts(struct('svg', true, 'png', false, 'fig_visible', C.fig_visible, 'size_px', C.fig_size_px), opts);

  files = [ dir(fullfile(dirs.raw, 'Curva_Polarizaci*_*.*DTA')) ];
  if isempty(files), return; end

  S = sort_tokens(files, 'polarization');
  xlsx = fullfile(dirs.processed, 'Polarization_results.xlsx');
  sheets = {};

  Pasc = struct('I_last_A',[],'V_last_V',[]);
  Pdsc = struct('I_last_A',[],'V_last_V',[]);
  meta_first = struct(); area_used = C.active_area_cm2_default; area_warn=false;

  fig_VI = figure('Name','POL V vs I'); set(fig_VI,'Visible',opts.fig_visible); hold on; grid on; xlabel('I (A)'); ylabel('V (V)');
  fig_VI_last = figure('Name','POL V vs I (last points)'); set(fig_VI_last,'Visible',opts.fig_visible); hold on; grid on; xlabel('I (A)'); ylabel('V (V)');
  fig_Vt = figure('Name','POL V vs t'); set(fig_Vt,'Visible',opts.fig_visible); hold on; grid on; xlabel('t (s)'); ylabel('V (V)');
  fig_It = figure('Name','POL I vs t'); set(fig_It,'Visible',opts.fig_visible); hold on; grid on; xlabel('t (s)'); ylabel('I (A)');
  fig_dVdI = figure('Name','POL dV/dI (processed)'); set(fig_dVdI,'Visible',opts.fig_visible); hold on; grid on; xlabel('I (A)'); ylabel('dV/dI (Ohm)');
  fig_delta = figure('Name','POL ΔV (DSC-ASC) vs I'); set(fig_delta,'Visible',opts.fig_visible); hold on; grid on; xlabel('I (A)'); ylabel('ΔV (V)');

  G = struct('global_idx', [], 't_s', [], 'V_V', [], 'I_A', []);
  gidx = 0;

  for k = 1:numel(S)
    fp = S(k).path; nm = S(k).name;
    Ctab = import_gamry_curves(fp);
    M = read_dta_metadata(fp);
    if k==1, meta_first = M; end
    if isfield(M,'area'), area_used = M.area; else, area_warn=true; end

    plot(fig_VI, Ctab.I_A, Ctab.V_V, '-', 'DisplayName', nm);
    plot(fig_Vt, Ctab.t_s, Ctab.V_V, '-', 'DisplayName', nm);
    plot(fig_It, Ctab.t_s, Ctab.I_A, '-', 'DisplayName', nm);

    cfg = C.step;
    P = last_points(Ctab.t_s, Ctab.I_A, Ctab.V_V, cfg);

    if S(k).dirflag
      Pasc.I_last_A = [Pasc.I_last_A; P.I_last_A];
      Pasc.V_last_V = [Pasc.V_last_V; P.V_last_V];
    else
      Pdsc.I_last_A = [Pdsc.I_last_A; P.I_last_A];
      Pdsc.V_last_V = [Pdsc.V_last_V; P.V_last_V];
    end

    n = numel(Ctab.t_s);
    G.global_idx = [G.global_idx; (gidx+1:gidx+n)'];
    G.t_s = [G.t_s; Ctab.t_s];
    G.V_V = [G.V_V; Ctab.V_V];
    G.I_A = [G.I_A; Ctab.I_A];
    gidx = gidx + n;
  end

  if area_warn
    warning('AREA not found in metadata for some files. Using default %.3g cm^2', area_used);
  end

  tol_abs = C.step.dsc_dup_tol_abs_A; tol_rel = C.step.dsc_dup_tol_rel;
  if ~isempty(Pdsc.I_last_A)
    mask = true(size(Pdsc.I_last_A));
    if numel(Pdsc.I_last_A) >= 2
      for m = 2:numel(Pdsc.I_last_A)
        a = Pdsc.I_last_A(m-1); b = Pdsc.I_last_A(m);
        if abs(b-a) <= max(tol_abs, tol_rel*max(1,abs(a)))
          mask(m) = false;
        end
      end
    end
    Pdsc.I_last_A = Pdsc.I_last_A(mask);
    Pdsc.V_last_V = Pdsc.V_last_V(mask);
  end

  plot(fig_VI_last, Pasc.I_last_A, Pasc.V_last_V, 'o-','DisplayName','ASC');
  plot(fig_VI_last, Pdsc.I_last_A, Pdsc.V_last_V, 's-','DisplayName','DSC'); legend(fig_VI_last,'show','Location','best');

  dVdI_asc = diff(Pasc.V_last_V) ./ max(1e-12, diff(Pasc.I_last_A));
  I_mid_asc = (Pasc.I_last_A(1:end-1) + Pasc.I_last_A(2:end))/2;
  dVdI_dsc = diff(Pdsc.V_last_V) ./ max(1e-12, diff(Pdsc.I_last_A));
  I_mid_dsc = (Pdsc.I_last_A(1:end-1) + Pdsc.I_last_A(2:end))/2;
  plot(fig_dVdI, I_mid_asc, dVdI_asc, 'o-','DisplayName','ASC');
  plot(fig_dVdI, I_mid_dsc, dVdI_dsc, 's-','DisplayName','DSC'); legend(fig_dVdI,'show','Location','best');

  if ~isempty(Pasc.I_last_A) && !isempty(Pdsc.I_last_A)
    Vd_interp = interp1(Pdsc.I_last_A, Pdsc.V_last_V, Pasc.I_last_A, 'linear', 'extrap');
    plot(fig_delta, Pasc.I_last_A, Vd_interp - Pasc.V_last_V, 'o-');
  end

  sheets(end+1,:) = {'META', {
    'DATE', safe_get(meta_first,'date',''); 'TIME', safe_get(meta_first,'time','');
    'ISTEP1', safe_get(meta_first,'istep1', NaN); 'TSTEP1', safe_get(meta_first,'tstep1', NaN);
    'ISTEP2', safe_get(meta_first,'istep2', NaN); 'TSTEP2', safe_get(meta_first,'tstep2', NaN);
    'SAMPLETIME', safe_get(meta_first,'sampletime', NaN);
    'AREA_cm2', area_used;
  }};
  sheets(end+1,:) = {'GLOBAL_raw', table_cells(G)};
  sheets(end+1,:) = {'ASC_processed', table_cells(struct('I_A',Pasc.I_last_A,'V_V',Pasc.V_last_V))};
  sheets(end+1,:) = {'DSC_processed', table_cells(struct('I_A',Pdsc.I_last_A,'V_V',Pdsc.V_last_V))};

  M = polarization_metrics(Pasc, Pdsc, area_used);
  sheets(end+1,:) = {'SUMMARY', {
    'n_ASC_pts', numel(Pasc.I_last_A);
    'n_DSC_pts', numel(Pdsc.I_last_A);
    'J_min_Acm2', M.J_range(1); 'J_max_Acm2', M.J_range(2);
    'V_min_V', M.V_range(1); 'V_max_V', M.V_range(2);
    'loop_area_V*A/cm2', M.loop_area_VAcm2;
  }};

  export_excel(xlsx, sheets);

  if opts.svg
    save_svg(fig_VI, fullfile(dirs.plots,'POL_VI.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_VI_last, fullfile(dirs.plots,'POL_VI_last.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_Vt, fullfile(dirs.plots,'POL_Vt.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_It, fullfile(dirs.plots,'POL_It.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_dVdI, fullfile(dirs.plots,'POL_dVdI.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_delta, fullfile(dirs.plots,'POL_deltaV.svg'), struct('png',opts.png,'size_px',opts.size_px));
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
    v = S.(f{j}); if iscell(v), v = v(:); else, v = num2cell(v(:)); end
    cells(2:end,j) = v;
  end
end


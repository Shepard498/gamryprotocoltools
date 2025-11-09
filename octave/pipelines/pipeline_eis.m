function pipeline_eis(dirs, opts)
  if nargin < 2, opts = struct(); end
  C = config();
  opts = merge_opts(struct('svg', true, 'png', false, 'fig_visible', C.fig_visible, 'size_px', C.fig_size_px), opts);

  files = dir(fullfile(dirs.raw, '*.*DTA'));
  if isempty(files), return; end

  stab = []; eis = [];
  for f = files'
    fp = fullfile(f.folder, f.name);
    if has_block(fp, 'ZCURVE')
      eis = [eis; f];
    elseif has_block(fp, 'CURVE') && ~isempty(regexp(f.name, 'EISPOT|EIS|Est|Stab', 'once','ignorecase'))
      stab = [stab; f];
    end
  end

  xlsx = fullfile(dirs.processed, 'EIS_results.xlsx');
  sheets = {};

  fig_stab = figure('Name','EIS Stabilization V vs t'); set(fig_stab,'Visible',opts.fig_visible); hold on; grid on; xlabel('t (s)'); ylabel('V (V)');
  for f = stab'
    fp = fullfile(f.folder, f.name);
    Ctab = import_gamry_curves(fp);
    plot(fig_stab, Ctab.t_s, Ctab.V_V, '-', 'DisplayName', f.name);
    sheets(end+1,:) = {sprintf('STAB_%s', sanitize_sheetname(f.name)), table_cells(Ctab)};
  end

  fig_nyq_3I = figure('Name','NYQ 1.5/9/18A'); set(fig_nyq_3I,'Visible',opts.fig_visible); hold on; grid on; axis equal; xlabel('Zre (Ohm)'); ylabel('-Zim (Ohm)');
  fig_nyq_0V = figure('Name','NYQ 0V (start/end)'); set(fig_nyq_0V,'Visible',opts.fig_visible); hold on; grid on; axis equal; xlabel('Zre (Ohm)'); ylabel('-Zim (Ohm)');

  fig_idc = figure('Name','EIS Idc vs index'); set(fig_idc,'Visible',opts.fig_visible); hold on; grid on; xlabel('Index'); ylabel('I_dc (A)');

  for f = eis'
    fp = fullfile(f.folder, f.name);
    Etab = import_gamry_eis(fp);
    M = read_dta_metadata(fp);

    fig_ind = figure('Name',['NYQ ' f.name]); set(fig_ind,'Visible',opts.fig_visible);
    plot(Etab.Zre_ohm, -Etab.Zim_ohm, '-o'); grid on; axis equal; xlabel('Zre (Ohm)'); ylabel('-Zim (Ohm)'); title(strrep(f.name,'_','\_'));
    if opts.svg, save_svg(fig_ind, fullfile(dirs.plots, [f.name '_nyquist.svg']), struct('png',opts.png,'size_px',opts.size_px)); end

    fig_bode = figure('Name',['Bode ' f.name]); set(fig_bode,'Visible',opts.fig_visible);
    subplot(2,1,1); loglog(Etab.freq_Hz, Etab.("|Z|_ohm"), '-'); grid on; ylabel('|Z| (Ohm)');
    subplot(2,1,2); semilogx(Etab.freq_Hz, Etab.phase_deg, '-'); grid on; xlabel('f (Hz)'); ylabel('Phase (deg)');
    if opts.svg, save_svg(fig_bode, fullfile(dirs.plots, [f.name '_bode.svg']), struct('png',opts.png,'size_px',opts.size_px)); end

    if is_zeroV(M)
      plot(fig_nyq_0V, Etab.Zre_ohm, -Etab.Zim_ohm, '-o','DisplayName', f.name);
    else
      plot(fig_nyq_3I, Etab.Zre_ohm, -Etab.Zim_ohm, '-o','DisplayName', label_current(M, Etab));
    end

    plot(fig_idc, 1:numel(Etab.Idc_A), Etab.Idc_A, '-','DisplayName', f.name);
    tgt = target_current_from_meta(M);
    if ~isnan(tgt), yline(fig_idc.CurrentAxes, tgt, '--'); end

    sheets(end+1,:) = {sprintf('EIS_%s', sanitize_sheetname(f.name)), table_cells(struct( ...
      'index', (1:numel(Etab.freq_Hz))', 'freq_Hz', Etab.freq_Hz, 'Zre_ohm', Etab.Zre_ohm, 'Zim_ohm', Etab.Zim_ohm, ...
      'Zphz_deg', Etab.phase_deg, 'Zmod_ohm', Etab.("|Z|_ohm"), 'Idc_A', Etab.Idc_A, 'Vdc_V', Etab.Vdc_V)) };
  end

  overview = {'Block','Name'};
  for f = stab', overview(end+1,:) = {'STAB', f.name}; end
  for f = eis',  overview(end+1,:) = {'EIS',  f.name}; end
  sheets = [{'META', overview}; sheets];

  export_excel(xlsx, sheets);

  if opts.svg
    save_svg(fig_stab, fullfile(dirs.plots,'EIS_stabilization_Vt.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_idc, fullfile(dirs.plots,'EIS_Idc_index.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_nyq_3I, fullfile(dirs.plots,'EIS_nyquist_3currents.svg'), struct('png',opts.png,'size_px',opts.size_px));
    save_svg(fig_nyq_0V, fullfile(dirs.plots,'EIS_nyquist_0V.svg'), struct('png',opts.png,'size_px',opts.size_px));
  end
end

function tf = has_block(fp, key)
  tf = false; fid = fopen(fp,'rt'); if fid<0, return; end; c = onCleanup(@() fclose(fid));
  while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    if strncmpi(strtrim(ln), key, numel(key)), tf = true; return; end
  end
end

function tf = is_zeroV(M)
  tf = isfield(M,'vdc') && abs(M.vdc) < 1e-6;
end

function s = label_current(M, E)
  if isfield(M,'i_dc'), s = sprintf('%.3g A', M.i_dc); return; end
  if any(~isnan(E.Idc_A)), s = sprintf('%.3g A', median(E.Idc_A(~isnan(E.Idc_A)))); return; end
  s = 'EIS';
end

function t = target_current_from_meta(M)
  t = NaN; if isfield(M,'i_dc'), t = M.i_dc; end
end

function cells = table_cells(S)
  f = fieldnames(S);
  hdr = f(:)'; n = numel(S.(f{1})); cells = [hdr; cell(n, numel(f))];
  for j = 1:numel(f)
    v = S.(f{j}); if iscell(v), v = v(:); else, v = num2cell(v(:)); end
    cells(2:end,j) = v;
  end
end


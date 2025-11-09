function pipeline_ocp(dirs, opts)
  if nargin < 2, opts = struct(); end
  C = config();
  opts = merge_opts(struct('svg', true, 'png', false, 'fig_visible', C.fig_visible, 'size_px', C.fig_size_px), opts);

  files = [ dir(fullfile(dirs.raw, 'OCP_*.*DTA')) ];
  if isempty(files), return; end

  xlsx = fullfile(dirs.processed, 'OCP_results.xlsx');
  sheets = {};

  for f = files'
    fp = fullfile(f.folder, f.name);
    Ctab = import_gamry_curves(fp);

    fig = figure('Name',['OCP ' f.name]); set(fig,'Visible',opts.fig_visible);
    plot(Ctab.t_s, Ctab.V_V, '-'); grid on; xlabel('t (s)'); ylabel('V (V)'); title(strrep(f.name,'_','\_'));
    if opts.svg
      save_svg(fig, fullfile(dirs.plots, [f.name '_ocp.svg']), struct('png',opts.png,'size_px',opts.size_px));
    end

    sheets(end+1,:) = {sprintf('OCP_%s', sanitize_sheetname(f.name)), table_cells(Ctab)};
  end

  export_excel(xlsx, sheets);
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


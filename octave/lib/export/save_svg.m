function save_svg(fig, outpath, opts)
  if nargin < 3, opts = struct(); end
  if ~isfield(opts,'png'), opts.png = false; end
  set(fig, 'PaperPositionMode','auto');
  print(fig, outpath, '-dsvg');
  if opts.png
    [p,n] = fileparts(outpath);
    print(fig, fullfile(p, [n '.png']), '-dpng', '-r300');
  end
end


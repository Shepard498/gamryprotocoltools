function save_svg(fig, outpath, opts)
  arguments
    fig (1,1)
    outpath (1,:) char
    opts.png (1,1) logical = false
    opts.size_px (1,2) double = [1600 900]
  end
  set(fig, 'Units','pixels', 'Position',[100 100 opts.size_px]);
  set(fig, 'PaperPositionMode','auto');
  print(fig, outpath, '-dsvg');
  if opts.png
    [p,n] = fileparts(outpath);
    print(fig, fullfile(p, [n '.png']), '-dpng', '-r300');
  end
end

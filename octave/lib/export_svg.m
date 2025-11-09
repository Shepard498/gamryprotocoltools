function export_svg(fig, outpath, varargin)
% export_svg(fig, outpath [, 'Width', 1600, 'Height', 1000])
% Saves the figure as vector SVG. Width/Height are in pixels (canvas size).
  p = inputParser;
  addParameter(p,'Width',1600,@(x)isnumeric(x)&&isscalar(x));
  addParameter(p,'Height',1000,@(x)isnumeric(x)&&isscalar(x));
  parse(p,varargin{:});
  [folder,~,~] = fileparts(outpath);
  if ~isempty(folder) && ~exist(folder,'dir'), mkdir(folder); end

  old = get(fig, 'paperpositionmode'); set(fig,'paperpositionmode','auto');
  try
    % -S sets canvas size; works for qt toolkit too.
    print(fig, outpath, '-dsvg', sprintf('-S%d,%d', p.Results.Width, p.Results.Height));
  finally
    set(fig,'paperpositionmode',old);
  end
end

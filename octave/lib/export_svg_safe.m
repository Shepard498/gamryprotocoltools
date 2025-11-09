function export_svg_safe(fig, outfn, varargin)
%EXPORT_SVG_SAFE  Robust SVG export with painters/gnuplot fallback.
%   Mirrors the pipeline-local helper previously bundled with the
%   polarización pipeline. Attempts a painters render first and falls back
%   to temporarily switching to the gnuplot toolkit when Octave/GL2PS
%   glitches occur. Errors are rethrown after restoring the original
%   toolkit where possible.

  try
    if ishghandle(fig)
      try, set(fig, 'renderer', 'painters'); catch, end
      drawnow();
    end
    export_svg_scaled(fig, outfn, varargin{:});
  catch
    oldtk = '';
    try
      oldtk = graphics_toolkit();
      graphics_toolkit('gnuplot');
      drawnow();
      export_svg_scaled(fig, outfn, varargin{:});
      if ~isempty(oldtk)
        graphics_toolkit(oldtk);
      end
    catch ME
      if ~isempty(oldtk)
        try, graphics_toolkit(oldtk); catch, end
      end
      rethrow(ME);
    end
  end
end

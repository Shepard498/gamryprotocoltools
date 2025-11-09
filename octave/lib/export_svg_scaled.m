function export_svg_scaled(hfig, filepath, varargin)
%EXPORT_SVG_SCALED  Clone a figure, scale visual elements, force a font, export SVG.
%
%   export_svg_scaled(h, 'out.svg', 'Scale',1.8, 'Width',1600, 'Height',1000,
%                     'LegendScale',1.8, 'LegendPad',[10 6], 'LegendCols',2,
%                     'LegendOutside',true)
%
% Notes
% - On-screen figure is untouched (we clone).
% - Typeface is enforced ONLY in the export (set FONT_FAMILY below).
% - Legend box gets extra space via LegendScale/LegendPad and optional columns.
% - Requires export_svg.m (wrapper) on path. Falls back to print -dsvg on error.
%
% Safer defaults (avoid spaces in font family for SVG)
FONT_FAMILY  = 'DejaVu Sans';
FONT_WEIGHT  = 'normal';
MIN_FONTSIZE = 4;

% ------------ params ------------
p = inputParser;
p.addParameter('Scale',          1.4,      @(x)isnumeric(x)&&x>0);
p.addParameter('Width',          1600,   @(x)isnumeric(x)&&x>0);
p.addParameter('Height',         1000,   @(x)isnumeric(x)&&x>0);
% legend controls
p.addParameter('LegendScale',    [],     @(x) isempty(x) || (isscalar(x)&&x>0));
p.addParameter('LegendPad',      [8 4],  @(x)isnumeric(x)&&numel(x)==2);   % [pxW pxH]
p.addParameter('LegendCols',     [],     @(x) isempty(x) || (isscalar(x)&&x>=1));
p.addParameter('LegendOutside',  false,  @islogical);
p.parse(varargin{:});
S = p.Results.Scale;
if isempty(p.Results.LegendScale), Ls = S; else, Ls = p.Results.LegendScale; end

% ------------ clone -------------
hclone = copyobj(hfig, 0);
set(hclone, 'Visible','off');

% ------------ scale lines / markers ------------
L = findobj(hclone, 'Type','line');
for k=1:numel(L)
  if isprop(L(k),'LineWidth'),  set(L(k),'LineWidth',  max(0.5, get(L(k),'LineWidth')  * S)); end
  if isprop(L(k),'MarkerSize'), set(L(k),'MarkerSize', max(1.0, get(L(k),'MarkerSize') * S)); end
end
G = findobj(hclone, '-regexp','Type','(scatter|hggroup|patch|bar)');
for k=1:numel(G)
  if isprop(G(k),'LineWidth'),  set(G(k),'LineWidth',  max(0.5, get(G(k),'LineWidth')  * S)); end
  if isprop(G(k),'MarkerSize'), set(G(k),'MarkerSize', max(1.0, get(G(k),'MarkerSize') * S)); end
end

% ------------ axes (ticks, box) ------------
AX = findobj(hclone, 'Type','axes');
for k=1:numel(AX)
  if isprop(AX(k),'FontSize'),  set(AX(k),'FontSize',  max(MIN_FONTSIZE, get(AX(k),'FontSize')*S)); end
  if isprop(AX(k),'LineWidth'), set(AX(k),'LineWidth', max(0.5, get(AX(k),'LineWidth')*S)); end
  if isprop(AX(k),'FontName'),  set(AX(k),'FontName',  FONT_FAMILY); end
  if isprop(AX(k),'FontWeight'),set(AX(k),'FontWeight',FONT_WEIGHT); end
  try
    enforce_text(get(AX(k),'XLabel'),S,FONT_FAMILY,FONT_WEIGHT,MIN_FONTSIZE);
    enforce_text(get(AX(k),'YLabel'),S,FONT_FAMILY,FONT_WEIGHT,MIN_FONTSIZE);
    enforce_text(get(AX(k),'ZLabel'),S,FONT_FAMILY,FONT_WEIGHT,MIN_FONTSIZE);
    enforce_text(get(AX(k),'Title'), S,FONT_FAMILY,FONT_WEIGHT,MIN_FONTSIZE);
  catch
  end
end

% ------------ legends (text + box) ------------
LG = findobj(hclone, 'Type','legend');
for k=1:numel(LG)
  % basic text/line scaling
  if isprop(LG(k),'Interpreter'), set(LG(k),'Interpreter','none'); end
  if isprop(LG(k),'FontSize'),    set(LG(k),'FontSize',   max(MIN_FONTSIZE, get(LG(k),'FontSize')*S)); end
  if isprop(LG(k),'LineWidth'),   set(LG(k),'LineWidth',  max(0.5, get(LG(k),'LineWidth')*S)); end
  if isprop(LG(k),'FontName'),    set(LG(k),'FontName',   FONT_FAMILY); end
  if isprop(LG(k),'FontWeight'),  set(LG(k),'FontWeight', FONT_WEIGHT); end

  % sanitize legend strings for SVG
  try
    strs = get(LG(k),'String');
    if iscell(strs), strs = cellfun(@svg_safe_text, strs, 'uni', 0); else, strs = svg_safe_text(strs); end
    set(LG(k),'String',strs);
  catch
  end

  % optional columns
  if ~isempty(p.Results.LegendCols) && isprop(LG(k),'NumColumns')
    set(LG(k),'NumColumns', p.Results.LegendCols);
  elseif isprop(LG(k),'NumColumns')
    % heuristic: if many entries, split in 2 columns
    try
      nlab = numel(get(LG(k),'String'));
      if nlab >= 8, set(LG(k),'NumColumns', 2); end
    catch
    end
  end

  % place outside if requested or if crowded
  if p.Results.LegendOutside && isprop(LG(k),'Location')
    set(LG(k),'Location','bestoutside');
  else
    % leave wherever it is (inside) but grow the bounding box
    try
      % 1) grow token size (MATLAB has ItemTokenSize; guard for Octave)
      if isprop(LG(k),'ItemTokenSize')
        ts = get(LG(k),'ItemTokenSize');
        if isnumeric(ts) && numel(ts)==2
          set(LG(k),'ItemTokenSize', max(1, round(ts*Ls)) );
        end
      end
      % 2) enlarge Position rectangle
      set(LG(k),'Units','pixels');
      pos = get(LG(k),'Position');  % [x y w h] in figure pixels (Octave honors pixels)
      pad = double(p.Results.LegendPad(:)');
      if numel(pad)~=2, pad = [8 4]; end
      pos(3) = pos(3)*Ls + pad(1);  % width
      pos(4) = pos(4)*Ls + pad(2);  % height
      set(LG(k),'Position',pos);
      % normalize back to avoid surprises downstream
      set(LG(k),'Units','normalized');
    catch
    end
  end
end

% ------------ all text (annotations, etc.) ------------
TX = findobj(hclone, 'Type','text');
for k=1:numel(TX)
  enforce_text(TX(k), S, FONT_FAMILY, FONT_WEIGHT, MIN_FONTSIZE);
end

% ------------ final safety tweaks ------------
set(hclone,'PaperPositionMode','auto');
try, set(hclone,'Renderer','painters'); end  % prefer painters for SVG

% ------------ export with fallback ------------
try
  export_svg(hclone, filepath, 'Width', p.Results.Width, 'Height', p.Results.Height);
catch
  % Fallback: direct print, which sometimes copes better with plotyy figures
  try
    print(hclone, filepath, '-dsvg', sprintf('-S%d,%d', p.Results.Width, p.Results.Height));
  catch ME2
    delete(hclone);
    rethrow(ME2);
  end
end

delete(hclone);
end

% ---- helper: enforce font + size on a text-like object ----
function enforce_text(h, S, family, weight, minfs)
  if isempty(h) || ~ishandle(h), return; end
  if isprop(h,'Interpreter'), set(h,'Interpreter','none'); end
  if isprop(h,'String')
    try, set(h,'String', svg_safe_text(get(h,'String'))); end
  end
  if isprop(h,'FontSize'),   set(h,'FontSize',   max(minfs, get(h,'FontSize')*S)); end
  if isprop(h,'FontName'),   set(h,'FontName',   family); end
  if isprop(h,'FontWeight'), set(h,'FontWeight', weight); end
end

% ---- helper: sanitize strings for SVG (escape & and remove odd chars) ----
function s = svg_safe_text(s)
  if iscell(s), s = cellfun(@svg_safe_text, s, 'uni', 0); return; end
  if ~ischar(s), return; end
  % Replace ampersands (XML) and normalize whitespace/newlines
  s = strrep(s, '&', 'and');
  s = regexprep(s, '[\r\n\t]+', ' ');
  % Optional: transliterate a few common accented letters to ASCII (gl2ps can be picky)
  repl = {'á','a';'é','e';'í','i';'ó','o';'ú','u';'Á','A';'É','E';'Í','I';'Ó','O';'Ú','U';'ñ','n';'Ñ','N'};
  for i=1:size(repl,1), s = strrep(s, repl{i,1}, repl{i,2}); end
end


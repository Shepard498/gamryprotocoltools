function [t, V, I] = import_gamry_dta(fp, opts)
%IMPORT_GAMRY_DTA  Robust .DTA importer for chrono-like files (t,V,I).
%   [t,V,I] = import_gamry_dta(filepath, opts)
% opts.header_guess_lines (default 160)
% opts.min_numeric_cols   (default 5)
% opts.decimal_comma_to_dot (default true)

  if nargin < 2, opts = struct(); end
  opts = merge_opts(struct('header_guess_lines',160, ...
                           'min_numeric_cols',  5,  ...
                           'decimal_comma_to_dot', true), opts);

  % Read entire file at once and normalize decimals
  try
    content = fileread(fp);
  catch
    error('No se pudo abrir: %s', fp);
  end
  
  if opts.decimal_comma_to_dot
    content = normalize_decimals(content);
  end
  
  % Split into lines and parse numeric data
  lines = regexp(content, '\r\n|\n|\r', 'split');
  data = zeros(0, opts.min_numeric_cols);
  
  for ln = lines'
    nums = sscanf(ln{1}, '%f');
    if numel(nums) >= opts.min_numeric_cols
      data(end+1,1:opts.min_numeric_cols) = nums(1:opts.min_numeric_cols)'; %#ok<AGROW>
    end
  end

  if isempty(data)
    error('No se encontraron líneas numéricas en: %s', fp);
  end
  if size(data,2) < 4
    error('El .DTA tiene menos de 4 columnas numéricas: %s', fp);
  end

  % Gamry (t,V,I) typical columns for chrono runs:
  t = data(:,2);
  V = data(:,3);
  I = data(:,4);
end


function str = normalize_decimals(str)
%NORMALIZE_DECIMALS Convert decimal commas to dots in a string or cell array
%   Centralizes decimal normalization logic used across import functions
  if isempty(str), return; end
  if iscell(str)
    str = cellfun(@normalize_decimals, str, 'UniformOutput', false);
  else
    str = strrep(str, ',', '.');
  end
end
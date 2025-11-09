function s = sanitize_sheetname(name)
  s = regexprep(name, '[/\\\?\*\[\]:]', '_');
  if numel(s) > 31, s = s(1:31); end
  if isempty(s), s = 'Sheet1'; end
end

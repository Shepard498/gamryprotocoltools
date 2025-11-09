function s = sanitize_sheet_name(s)
%SANITIZE_SHEET_NAME  Excel-safe sheet name: <=31 chars, invalid chars replaced.
  if isempty(s), s = 'Sheet'; end
  bad = '[]:*?/\"';
  for kk = 1:numel(bad)
    s = strrep(s, bad(kk), '_');
  end
  s = strtrim(s);
  if numel(s) > 31
    s = s(1:31);
  end
  if isempty(s)
    s = 'Sheet';
  end
end


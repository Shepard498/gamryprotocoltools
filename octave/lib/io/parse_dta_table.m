function T = parse_dta_table(filepath, key)
  fid = fopen(filepath, 'rt');
  if fid < 0, error('Cannot open %s', filepath); end
  cleaner = onCleanup(@() fclose(fid));

  headerline = '';
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if strncmpi(strtrim(ln), key, numel(key))
      headerline = ln; break;
    end
  end
  if isempty(headerline), error('Key %s not found in %s', key, filepath); end

  hdr = '';
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if ~isempty(strtrim(ln))
      hdr = ln; break;
    end
  end
  if isempty(hdr), error('No column header after %s in %s', key, filepath); end

  headers = regexp(strtrim(hdr), '[\t,;\s]+', 'split');

  data = [];
  while true
    pos = ftell(fid);
    ln = fgetl(fid);
    if ~ischar(ln) || isempty(strtrim(ln)), break; end
    toks = regexp(strtrim(ln), '[\t,;\s]+', 'split');
    vals = str2num_locale(toks);
    if any(isnan(vals))
      fseek(fid, pos, 'bof');
      break;
    end
    data(end+1,1:numel(vals)) = vals; %#ok<AGROW>
  end

  headers = lower(strtrim(headers));
  headers = regexprep(headers, '[^a-z0-9_]+', '_');
  T = struct();
  for j = 1:numel(headers)
    col = data(:, j);
    T.(headers{j}) = col;
  end
end


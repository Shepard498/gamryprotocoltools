function T = parse_dta_table(filepath, key)
  % Robust Gamry DTA parser that tolerates decimal commas.
  fid = fopen(filepath, 'rt');
  if fid < 0, error('Cannot open %s', filepath); end
  cleaner = onCleanup(@() fclose(fid));

  % 1) Seek the section header line starting with KEY (e.g., CURVE or ZCURVE)
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if strncmpi(strtrim(ln), key, numel(key)), break; end
  end
  if ~ischar(ln), error('Key %s not found in %s', key, filepath); end

  % 2) Next non-empty line is the column header
  hdr = '';
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if ~isempty(strtrim(ln)), hdr = ln; break; end
  end
  if isempty(hdr), error('No column header after %s in %s', key, filepath); end

  headers = regexp(strtrim(hdr), '[\t,; \s]+', 'split');
  ncols   = numel(headers);
  data    = [];

  % 3) Read numeric rows; tolerate commas and ragged columns. Stop on blank or next section.
  while true
    pos = ftell(fid);
    ln  = fgetl(fid);
    if ~ischar(ln), break; end
    s = strtrim(ln);
    if isempty(s), break; end
    % If we hit another section header (ALLCAPS word), stop.
    if ~isempty(regexp(s, '^[A-Z][A-Z0-9_]+', 'once')) && ~isempty(strfind(s, ':'))
      fseek(fid, pos, 'bof'); break;
    end
    toks = regexp(s, '[\t,; \s]+', 'split');
    vals = str2num_locale(toks);
    if all(isnan(vals))  % non-numeric fluff line, skip
      continue;
    end
    % Normalize width: pad/truncate to header count
    if numel(vals) < ncols, vals(end+1:ncols) = NaN; end
    if numel(vals) > ncols, vals = vals(1:ncols); end
    data(end+1,1:ncols) = vals; %#ok<AGROW>
  end

  if isempty(data)
    error('No numeric rows after %s in %s (check decimal commas & headers).', key, filepath);
  end

  % 4) Normalize header names and build struct
  headers = lower(strtrim(headers));
  headers = regexprep(headers, '[^a-z0-9_]+', '_');
  T = struct();
  for j = 1:ncols
    T.(headers{j}) = data(:, j);
  end
end


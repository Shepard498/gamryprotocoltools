function T = parse_dta_table(filepath, key)
  % PARSE_DTA_TABLE Generic text parser for Gamry .DTA tables.
  % key: 'CURVE' (chrono/polarization), 'ZCURVE' (EIS), etc.
  % Returns struct with fields = lowercase normalized column names.
  fid = fopen(filepath, 'rt');
  if fid < 0, error('Cannot open %s', filepath); end
  cleaner = onCleanup(@() fclose(fid));


  % Find the header line that starts with key
  headerline = '';
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if strncmpi(strtrim(ln), key, numel(key))
      headerline = ln; break;
    end
  end
  if isempty(headerline), error('Key %s not found in %s', key, filepath); end


  % The next non-empty line is assumed to be column headers
  hdr = '';
  while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    if ~isempty(strtrim(ln))
      hdr = ln; break;
    end
  end
  if isempty(hdr), error('No column header after %s in %s', key, filepath); end


  % Split headers by common delimiters
  delims = {"\t", ",", ";", "\s+"};
  for d = 1:numel(delims)
    h = regexp(hdr, delims{d}, 'split');
    if numel(h) > 1, delim = delims{d}; headers = h; break; end
  end
  if ~exist('headers','var'), headers = {hdr}; delim = "\s+"; end


  % Read numeric block until a blank or non-numeric line appears
  data = [];
  while true
    pos = ftell(fid);
    ln = fgetl(fid);
    if ~ischar(ln) || isempty(strtrim(ln))
      break;
    end
    nums = regexp(strtrim(ln), delim, 'split');
    vals = str2double(nums);
    if any(isnan(vals))
      fseek(fid, pos, 'bof');
    break;
    end
    data(end+1,1:numel(vals)) = vals; %#ok<AGROW>
  end


  % Build struct with normalized field names
  headers = lower(strtrim(headers));
  headers = regexprep(headers, '[^a-z0-9_]+', '_');
  T = struct();
  for j = 1:numel(headers)
    col = data(:, j);
    T.(headers{j}) = col;
  end
end

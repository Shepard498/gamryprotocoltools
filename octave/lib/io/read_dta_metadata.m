function M = read_dta_metadata(filepath)
  M = struct();
  fid = fopen(filepath, 'rt');
  if fid < 0, error('Cannot open %s', filepath); end
  cleaner = onCleanup(@() fclose(fid));
  while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    if contains(ln, '\t')
      toks = regexp(ln, '\t', 'split');
      if numel(toks) >= 3
        key = strtrim(toks[1]); typ = strtrim(toks{2}); val = strtrim(toks{3}); %#ok<NASGU>
        val = regexprep(val, '(?<=[0-9]),(?=[0-9])', '.');
        numv = str2double(val);
        if ~isnan(numv)
          M.(sanitize_key(key)) = numv;
        else
          M.(sanitize_key(key)) = val;
        end
      end
    end
  end
end

function k = sanitize_key(s)
  k = lower(regexprep(s, '[^a-z0-9_]+', '_'));
end


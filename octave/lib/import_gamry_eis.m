function [f, Zre, Zim, V, I] = import_gamry_eis(filepath, imp)
% Header-aware EIS import for Gamry ZCURVE tables (comma-safe paths).
  if nargin<2 || ~isstruct(imp), imp = struct(); end
  if ~isfield(imp,'decimal_comma_to_dot'), imp.decimal_comma_to_dot = true; end

  assert(exist(filepath,'file')==2, 'Cannot find file: %s', filepath);

  % --- read whole file as text (don’t touch the path!) ---
  txt = fileread(filepath);                 % avoids fopen issues with commas
  if imp.decimal_comma_to_dot
    txt = strrep(txt, ',', '.');            % normalize ONLY the contents
  end

  % --- split into lines ---
  C = regexp(txt, '\r\n|\n|\r', 'split');

  % --- find ZCURVE header & names ---
  iZ = find(~cellfun('isempty', regexp(C,'^ZCURVE')), 1, 'first');
  assert(~isempty(iZ), 'ZCURVE block not found in %s', filepath);

  hdrLine = strtrim(C{iZ+1});
  toks = regexp(hdrLine, '\s+', 'split');
  if isempty(toks{1}), toks(1)=[]; end

  idxMap = containers.Map('KeyType','char','ValueType','double');
  for j=1:numel(toks), idxMap(toks{j}) = j; end

  need = {'Freq','Zreal','Zimag'};
  for j=1:numel(need)
    assert(isKey(idxMap, need{j}), 'Column "%s" not found in %s', need{j}, filepath);
  end
  hasIdc = isKey(idxMap,'Idc');  hasVdc = isKey(idxMap,'Vdc');

  % --- parse numeric rows ---
  D = [];
  for k = iZ+3 : numel(C)
    line = strtrim(C{k});
    if isempty(line), break; end
    nums = sscanf(line, '%f');
    if isempty(nums), break; end
    if numel(nums) >= numel(toks)
      row = nums(end-numel(toks)+1:end).';
      D(end+1,1:numel(row)) = row; %#ok<AGROW>
    end
  end
  assert(~isempty(D), 'No numeric ZCURVE rows parsed in %s', filepath);

  % --- extract vectors ---
  f   = D(:, idxMap('Freq'));
  Zre = D(:, idxMap('Zreal'));
  Zim = D(:, idxMap('Zimag'));
  if hasIdc, I = D(:, idxMap('Idc')); else, I = NaN(size(f)); end
  if hasVdc, V = D(:, idxMap('Vdc')); else, V = NaN(size(f)); end
end


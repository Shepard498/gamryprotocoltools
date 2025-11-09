function files = list_files(folder, pattern)
%LIST_FILES  List *.DTA in folder and parse tokens using regex pattern.
%   files is struct array with:
%     .name, .full, .tokens (cellstr of captures)
%   Example pattern for Polarización:
%     '^Curva_Polarizacion_(Asc|Dsc)_(.+?)_(.+?)_#(\d+)\.DTA$'

  if nargin < 2 || isempty(pattern)
    pattern = '^(.*)$'; % match anything
  end
  L = dir(fullfile(folder, '*.DTA'));
  files = struct('name',{},'full',{},'tokens',{});
  for kk = 1:numel(L)
    fn = L(kk).name;
    tok = regexp(fn, pattern, 'tokens', 'once');
    if isempty(tok), continue; end
    files(end+1).name   = fn; %#ok<AGROW>
    files(end).full     = fullfile(folder, fn);
    files(end).tokens   = tok;
  end
end


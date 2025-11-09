function [hdr, cells] = struct2cellsheet(S)
  f = fieldnames(S);
  hdr = f(:)';
  n = numel(S.(f{1}));
  cells = cell(n, numel(f));
  for j = 1:numel(f)
    v = S.(f{j});
    if iscell(v), v = v(:); else, v = num2cell(v(:)); end
    cells(:,j) = v;
  end
end

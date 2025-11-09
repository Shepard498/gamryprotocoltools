function cells = struct2cellsheet(S)
  f = fieldnames(S);
  hdr = f(:)';
  n = numel(S.(f{1}));
  cells = [hdr; cell(n, numel(f))];
  for j = 1:numel(f)
    v = S.(f{j});
    if iscell(v), v = v(:); else, v = num2cell(v(:)); end
    cells(2:end,j) = v;
  end
end


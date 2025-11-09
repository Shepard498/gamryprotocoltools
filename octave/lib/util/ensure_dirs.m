function ensure_dirs(dirs)
  f = fieldnames(dirs);
  for k = 1:numel(f)
    d = dirs.(f{k});
    if ~exist(d, 'dir'), mkdir(d); end
  end
end

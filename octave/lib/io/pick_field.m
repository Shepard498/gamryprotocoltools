function v = pick_field(T, fn, candidates)
  v = [];
  for k = 1:numel(candidates)
    name = lower(regexprep(candidates{k}, '[^a-z0-9_]+','_'));
    h = fn(strcmpi(fn, name));
    if ~isempty(h), v = T.(h{1}); return; end
  end
end

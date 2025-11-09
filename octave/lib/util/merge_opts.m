function out = merge_opts(def, opt)
  if nargin < 2 || isempty(opt), out = def; return; end
  out = def;
  k = fieldnames(opt);
  for i = 1:numel(k)
    out.(k{i}) = opt.(k{i});
  end
end

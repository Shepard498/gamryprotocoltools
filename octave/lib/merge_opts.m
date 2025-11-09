function out = merge_opts(defaults, opts)
%MERGE_OPTS  Shallow-merge struct fields of opts into defaults.
%   out = MERGE_OPTS(defaults, opts)
% If a field exists in opts, it overwrites defaults. Non-destructive.
  if nargin < 2 || isempty(opts)
    out = defaults; return;
  end
  out = defaults;
  f = fieldnames(opts);
  for k = 1:numel(f)
    out.(f{k}) = opts.(f{k});
  end
end

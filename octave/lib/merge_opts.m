function out = merge_opts(base, extra)
%MERGE_OPTS  Shallow-merge two structs (extra overrides base).
  if nargin < 1 || isempty(base),  base  = struct(); end
  if nargin < 2 || isempty(extra), extra = struct(); end
  out = base;
  if ~isstruct(extra), error('extra must be a struct'); end
  f = fieldnames(extra);
  for kk = 1:numel(f)
    out.(f{kk}) = extra.(f{kk});
  end
end


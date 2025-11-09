function segments = detect_plateaus(I, N, opts)
%DETECT_PLATEAUS  Find N plateaus (sub-steps) in current trace I.
%   segments = [s1 e1; s2 e2; ... sN eN]
% opts fields:
%   smooth_frac   (default 0.01)  % moving-average window fraction of N
%   edge_guard    (default 0.05)  % ignore first/last fraction when finding splits
%   min_len_frac  (default 0.05)  % minimal segment length as fraction of N

  if nargin < 2 || isempty(N), N = 2; end
  if nargin < 3, opts = struct(); end
  L = numel(I);
  if L < max(10, 2*N), segments = [1, L]; return; end

  opts = merge_opts(struct('smooth_frac', 0.01, ...
                           'edge_guard',  0.05, ...
                           'min_len_frac',0.05), opts);

  w = max(3, ceil(opts.smooth_frac * L));
  if mod(w,2) == 0, w = w + 1; end

  % moving average (simple FIR)
  Ipad = [I(1)*ones(w,1); I(:); I(end)*ones(w,1)];
  If   = filter(ones(w,1)/w, 1, Ipad);
  If   = If(w+1:w+L);

  if N == 1
    segments = [1, L];
    return;
  end

  % find (N-1) split points by largest |diff|
  g = abs(diff(If));
  guard = max(1, round(opts.edge_guard * L));
  g(1:guard) = 0; g(end-guard+1:end) = 0;

  % pick top (N-1) peaks, enforce spacing
  [~, idx_sorted] = sort(g(:), 'descend');
  splits = [];
  for kk = 1:numel(idx_sorted)
    cand = idx_sorted(kk);
    ok = true;
    for jj = 1:numel(splits)
      if abs(cand - splits(jj)) < round(opts.min_len_frac * L)
        ok = false; break;
      end
    end
    if ok
      splits(end+1) = cand; %#ok<AGROW>
      if numel(splits) == (N-1), break; end
    end
  end

  if isempty(splits)
    % fallback: split at median index if nothing stands out
    splits = round(linspace(1, L-1, N+1));
    splits = splits(2:end-1);
  else
    splits = sort(splits);
  end

  % build segments from splits
  b = [1, splits+1];
  e = [splits, L];
  segments = [b(:), e(:)];
end


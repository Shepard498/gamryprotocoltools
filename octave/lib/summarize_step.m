function S = summarize_step(t, V, I, seg, mode, params)
%SUMMARIZE_STEP  Summarize a sub-step by last point or tail average.
%   S has fields: .I, .V, .t, .idx_mask
% mode: 'last' | 'tail_seconds' | 'tail_frac'
% params.tail_seconds (for 'tail_seconds'), params.tail_frac (for 'tail_frac')

  if nargin < 5 || isempty(mode), mode = 'last'; end
  if nargin < 6, params = struct(); end

  idx = seg(1):seg(2);
  tl  = t(idx); Vl = V(idx); Il = I(idx);

  switch lower(mode)
    case 'last'
      S.I = Il(end); S.V = Vl(end); S.t = tl(end);
      m = false(size(I)); m(idx(end)) = true; S.idx_mask = m;

    case 'tail_seconds'
      if ~isfield(params, 'tail_seconds') || isempty(params.tail_seconds)
        params.tail_seconds = 10; % default
      end
      t0 = max(tl) - params.tail_seconds;
      wid = (tl >= t0);
      if ~any(wid), wid(end) = true; end
      S.I = mean(Il(wid)); S.V = mean(Vl(wid)); S.t = mean(tl(wid));
      m = false(size(I)); m(idx(1)-1 + find(wid)) = true; S.idx_mask = m;

    case 'tail_frac'
      if ~isfield(params, 'tail_frac') || isempty(params.tail_frac)
        params.tail_frac = 0.30;
      end
      n  = numel(tl);
      n0 = max(1, floor((1-params.tail_frac) * n));
      wid = false(size(tl)); wid(n0:end) = true;
      S.I = mean(Il(wid)); S.V = mean(Vl(wid)); S.t = mean(tl(wid));
      m = false(size(I)); m(idx(1)-1 + find(wid)) = true; S.idx_mask = m;

    otherwise
      error('Unknown mode: %s', mode);
  end
end


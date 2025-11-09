function idx = step_boundaries(t, i, cfg)
  di = [0; abs(diff(i))./max(1e-9, diff(t))];
  thr = cfg.deriv_thr_frac * max(di);
  raw = [1; 1+find(di(2:end) > thr); numel(t)];
  idx = unique(raw);
end

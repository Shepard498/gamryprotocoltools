function P = last_points(t, i, v, cfg)
  idx = step_boundaries(t, i, cfg);
  I_last = []; V_last = []; T_last = []; step_idx = [];
  for k = 1:numel(idx)-1
    a = idx(k); b = idx(k+1)-1; if b<=a, continue; end
    if t(b)-t(a) < cfg.min_step_s, continue; end
    I_last(end+1,1) = i(b);
    V_last(end+1,1) = v(b);
    T_last(end+1,1) = t(b);
    step_idx(end+1,1) = k;
  end
  P = struct('I_last_A', I_last, 'V_last_V', V_last, 't_last_s', T_last, 'step_idx', step_idx);
end

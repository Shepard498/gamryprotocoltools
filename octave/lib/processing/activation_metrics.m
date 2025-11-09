function S = activation_metrics(global_last, area_cm2)
  if nargin < 2 || isempty(area_cm2), area_cm2 = 25; end
  S = struct();
  if isempty(global_last.I_A)
    S.n_steps = 0; S.n_cycles = 0; S.t_total_s = 0; S.I_range = [NaN NaN]; S.V_range = [NaN NaN]; return;
  end
  S.n_steps = numel(global_last.I_A);
  S.n_cycles = numel(unique(global_last.cycle));
  S.t_total_s = max(global_last.t_s) - min(global_last.t_s);
  S.I_range = [min(global_last.I_A) max(global_last.I_A)];
  S.V_range = [min(global_last.V_V) max(global_last.V_V)];
  S.J_range = S.I_range / area_cm2;
end

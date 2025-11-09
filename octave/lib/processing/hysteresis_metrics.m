function H = hysteresis_metrics(Pasc, Pdsc)
  if isempty(Pasc.I_last_A) || isempty(Pdsc.I_last_A)
    H = struct('mean_dV', NaN, 'max_dV', NaN); return;
  end
  Vd = interp1(Pdsc.I_last_A, Pdsc.V_last_V, Pasc.I_last_A, 'linear', 'extrap');
  dV = Vd - Pasc.V_last_V;
  H = struct('mean_dV', mean(dV), 'max_dV', max(abs(dV)));
end

function out = eis_metrics(E)
  N = min(10, numel(E.Zre_ohm));
  out = struct();
  out.R_ohm_est = min(E.Zre_ohm(1:N));
  out.freq_top_Hz = E.freq_Hz(1);
end

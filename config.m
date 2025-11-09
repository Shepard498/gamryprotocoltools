function C = config()
  % Global configuration (edit as needed)
  C = struct();
  % Electrode area [cm^2]. Will be overwritten per-file if found in DTA metadata.
  C.active_area_cm2_default = 25;  % warn if metadata missing

  % Figure defaults
  C.fig_size_px = [1600 900];
  C.fig_visible = 'on';  % 'off' for headless

  % Parsing
  C.decimal_policy = 'auto'; % 'auto' flips comma to dot when needed

  % Step detection & duplicates
  C.step = struct('min_step_s', 15, ...
                  'window_frac', 0.15, ...   % for averaging windows when needed
                  'deriv_thr_frac', 0.2, ... % threshold as fraction of max |di/dt|
                  'dsc_dup_tol_abs_A', 0.02, ...
                  'dsc_dup_tol_rel', 0.01);
end

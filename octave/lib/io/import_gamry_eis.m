function R = import_gamry_eis(filepath)
  T = parse_dta_table(filepath, 'ZCURVE');
  fn = fieldnames(T);
  get = @(names) pick_field(T, fn, names);

  f  = get({'freq_hz','frequency_hz','freq'});
  zr = get({'zreal_ohm','zre_ohm','z_real','zre'});
  zi = get({'zimag_ohm','zim_ohm','z_imag','zim'});
  zm = get({'|z|_ohm','zmod_ohm','z_magnitude'});
  ph = get({'phase_deg','phz_deg','z_phase_deg'});
  vdc= get({'vdc_v','vdcoffset_v','v_dc'});
  idc= get({'idc_a','idcoffset_a','i_dc'});

  if isempty(f) || isempty(zr) || isempty(zi)
    error('Unrecognized ZCURVE columns in %s', filepath);
  end
  if isempty(zm), zm = sqrt(zr.^2 + zi.^2); end
  if isempty(ph), ph = atan2d(zi, zr); end

  R = struct('freq_Hz', f(:), 'Zre_ohm', zr(:), 'Zim_ohm', zi(:), ...
             '|Z|_ohm', zm(:), 'phase_deg', ph(:), 'Vdc_V', vdc(:), 'Idc_A', idc(:));
end


function M = polarization_metrics(Pasc, Pdsc, area_cm2)
  if nargin < 3 || isempty(area_cm2), area_cm2 = 25; end
  Jasc = Pasc.I_last_A / area_cm2; Jdsc = Pdsc.I_last_A / area_cm2;
  if ~isempty(Jasc) && ~isempty(Jdsc)
    Vdsc_interp = interp1(Jdsc, Pdsc.V_last_V, Jasc, 'linear', 'extrap');
    dV = Vdsc_interp - Pasc.V_last_V;
    loop_area = trapz(Jasc, abs(dV));
  else
    dV = []; loop_area = NaN;
  end
  M = struct();
  M.Jasc_Acm2 = Jasc; M.Vasc_V = Pasc.V_last_V;
  M.Jdsc_Acm2 = Jdsc; M.Vdsc_V = Pdsc.V_last_V;
  M.dV_V = dV; M.loop_area_VAcm2 = loop_area;
  M.J_range = [min([Jasc;Jdsc]), max([Jasc;Jdsc])];
  M.V_range = [min([Pasc.V_last_V; Pdsc.V_last_V]), max([Pasc.V_last_V; Pdsc.V_last_V])];
end


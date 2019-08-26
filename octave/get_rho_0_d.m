#input: d dimension
#output: singlet
function rho_0_d = get_rho_0_d(d)
  rho_0_d = zeros(d);
  rho_0_d(1,1) = 1;
  rho_0_d(1,d) = 1;
  rho_0_d(d,1) = 1;
  rho_0_d(d,d) = 1;
  rho_0_d = rho_0_d / trace(rho_0_d);
endfunction
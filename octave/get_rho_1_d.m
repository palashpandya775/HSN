#input: d dimension
#output: identity noise
function rho_1_d = get_rho_1_d(d)
  rho_1_d = eye(d);
  rho_1_d = rho_1_d / trace(rho_1_d);
endfunction
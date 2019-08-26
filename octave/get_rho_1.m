#fixed d = 4
#input: none
#output: identity noise
function rho_1 = get_rho_1()
  rho_1 = eye(4);
  rho_1 = rho_1 / trace(rho_1);
endfunction
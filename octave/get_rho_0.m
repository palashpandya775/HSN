#fixed d = 4
#input: none
#output: singlet
function rho_0 = get_rho_0()
  rho_0 = zeros(4);
  rho_0(1,1) = 1;
  rho_0(1,4) = 1;
  rho_0(4,1) = 1;
  rho_0(4,4) = 1;
  rho_0 = rho_0 / trace(rho_0);
endfunction
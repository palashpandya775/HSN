function d = bures(rho_1, rho_2)
  F = (trace(mat_sqrt(mat_sqrt(rho_1)*rho_2*mat_sqrt(rho_1))))^2;
  d = 2*(1-sqrt(F));
endfunction
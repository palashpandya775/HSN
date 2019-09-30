#input: rho_0_d state of dimension d to test
#input: rho_1_d sep state of dimension d to start testing
#input: ds vector containing subsystems dimensions
#output: css_d closest seperable state of dimension d found
#assume: d = prod(ds)
function css_d = gilbert_d(rho_0_d, rho_1_d, ds)
  tic;
  CS_MAX = 5000;
  D_vals = [];
  
  #dimensions
  d = rows(rho_0_d);
  
  if(d != prod(ds))
    disp("error: dimensions mismatch");
    return;
  endif

  #disp("step 1");
  cs = 0;
  ct = 0;
  
  while(cs < CS_MAX)
    rho_2_d = 1;
    for i = ds
      A = randn(i,1) + randn(i,1)*i;
      A = A / norm(A);
      A = A * A';
      rho_2_d = kron(rho_2_d, A);
    endfor
    
    previous_d = trace((rho_0_d - rho_1_d)^2);
    
    #disp("step 2");
    pre_crit = trace((rho_2_d - rho_1_d)*(rho_0_d - rho_1_d));
    if(pre_crit <= 0)
      continue;
    endif
    
    #disp("step 3");
    #this step is skipped in this basic implementation
        
    #disp("step 4");
    sub_rho_2_rho_1 = rho_2_d-rho_1_d;
    min_p = -trace((sub_rho_2_rho_1)*(rho_0_d-rho_2_d))/trace((sub_rho_2_rho_1)*(sub_rho_2_rho_1));
    
    #disp("step 5");
    p = min_p;
    rho_1_candidate_d = p*rho_1_d + (1-p)*rho_2_d;
    result_d = trace((rho_0_d - rho_1_candidate_d)^2);
    if(result_d < previous_d && p >= 0 && p <=1)
      rho_1_d = rho_1_candidate_d;
      D_vals(end+1) = trace((rho_0_d - rho_1_d)^2);
      cs = cs + 1;
      
      printf("%d\t%f\n", cs, D_vals(end));
      fflush(stdout);
    endif
  
  #disp("step 6");
  endwhile
  css_d = rho_1_d;
  
  printf("ct: %d\n", ct);
  printf("cs: %d\n", cs);
  printf("\ntrace(css) = %f\n", trace(css_d));
  printf("eig(css) = %f\n", transpose(eig(css_d)));
  printf("\n");
  toc
endfunction

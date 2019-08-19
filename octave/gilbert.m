#Bianka is acknowledged hereby by Omer for reminding us all about the sqrt().
#Thank you Bianka.

function css = gilbert(rho_0, rho_1)
  h = 0.001;
  N = 1000;
  D_vals = [];
  
  #step 1
  cs = 0;
  ct = 0;
  while(ct < N)
    ct = ct + 1;
    
    A = rand(2,1);
    A = A / norm(A);
    
    B = rand(2,1);
    B = B / norm(B);
    
    A = A * A';
    B = B * B';
    
    rho_2 = kron(A,B);
    
    previous_d = trace((rho_0 - rho_1)^2);
    
    #disp("step 2");
    pre_crit = trace((rho_2 - rho_1)*(rho_0 - rho_1));
    if(pre_crit <= 0)
      continue;
    endif
    
    #disp("step 3");
    #this step is skipped in this basic implementation
    
    #disp("step 4");
    p = 0;
    min_p = trace((rho_0 - p*rho_1 - (1-p)*rho_2)^2);
    while(p <= 1)
      found_p = trace((rho_0 - p*rho_1 - (1-p)*rho_2)^2);
      if(found_p < min_p)
        min_p = found_p;
      endif
      p = p + h;
    endwhile
  
    #disp("step 5");
    p = min_p;
    rho_1_candidate = p*rho_1 + (1-p)*rho_2;
    result_d = trace((rho_0 - rho_1_candidate)^2);
    if(result_d < previous_d)
      rho_1 = rho_1_candidate;
      D_vals(end+1) = trace((rho_0 - rho_1)^2);
      cs = cs + 1
    endif
    
  #disp("step 6");
  endwhile
  
  disp(D_vals');
  css = rho_1;
endfunction
  
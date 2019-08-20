#Bianka is acknowledged hereby by Omer for reminding us all about the sqrt().
#Thank you Bianka.

function css = gilbert(rho_0, rho_1)
  tic;
  CT_MAX = 1000;
  h = 0.001;
  D_vals = [];
  
  #disp("step 1");
  cs = 0;
  ct = 0;
  while(ct < CT_MAX)
    ct = ct + 1;
    
    A = zeros(2,1);
    B = zeros(2,1);
    for k = [1:2]
      A(k) = randn() + ((-1)^(randi(2)-1)) * (randn()*i);
      B(k) = randn() + ((-1)^(randi(2)-1)) * (randn()*i);
    endfor
    A = A / norm(A);
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
      cs = cs + 1;
      
      printf("%d\t%f\n", cs, D_vals(end));
      fflush(stdout);
    endif
    
  #disp("step 6");
  endwhile
  css = rho_1;
  
  printf("\ntrace(css) = %f\n", trace(css));
  printf("eig(css) = %f\n", transpose(eig(css)));
  printf("\n");
  toc
endfunction

function rho_0 = get_rho_0()
  rho_0 = zeros(4);
  rho_0(1,1) = 1;
  rho_0(1,4) = 1;
  rho_0(4,1) = 1;
  rho_0(4,4) = 1;
  rho_0 = rho_0 / trace(rho_0);
endfunction

function rho_1 = get_rho_1()
  rho_1 = eye(4);
  rho_1 = rho_1 / trace(rho_1);
endfunction

#Return <i| vector 
#(indexed from 0)
function r = bra(i, d)
  r = zeros(1, d);
  r(i+1) = 1;
endfunction
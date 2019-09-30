#Return |i> vector
#(indexed from 0)
function r = ket(i, d)
  r = zeros(d, 1);
  r(i+1) = 1;
endfunction
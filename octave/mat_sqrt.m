#My matrix square root with exception handling when argument of sqrt is smaller than epsilon
function msqrt = mat_sqrt(M)
  [eigenvectors, eigenvalues] = eig(M);
  msqrt = zeros(size(M));
  for i = 1:columns(eigenvalues)
    if(eigenvalues(i,i) <= eps)
      temp_sqrt = 0;
    else
      temp_sqrt = sqrt(eigenvalues(i,i));
      msqrt = msqrt + temp_sqrt * eigenvectors(:,i) * eigenvectors(:,i)';
    endif
  endfor
endfunction
function res = reflex(original)

  N = length(original);
  assert(N>2, 'signal too short')

  res(1) = original(1);
  res(N) = abs(original(N));
  for i=2:N-1
    res(i) = original(i);
    res(2*N-i) = conj(original(i)); % Si se suma 1 al argumento de res, se duplica el último valor
  end
  assert(abs(original(N))==res(N), ...
         'original block inside res is not complete')
  assert(real(res(N))==res(N), ...
         'Middle element must be real')
end

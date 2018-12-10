function [ejet,p,Fs] = myifft(Y_reflected, df)

  L = length(Y_reflected);

  assert(mod(L,2)==0, 'Only even lengths, by now')
  N = (L/2)+1;

  try
    assert(isequal(fliplr(Y_reflected(N+1:L)), ...
                   conj(Y_reflected(2:N-1))), ...
           'Y_reflected problem')
  catch
    assert(check_reflection_onebyone(Y_reflected), ...
           'Y_reflected problem')
  end

  p = ifft(Y_reflected);
  assert(mean(abs(imag(p))) ...
         < mean(abs(real(p)))/100, ...
         'Nonzero imaginary part')
  p = real(p);
  
  % Normalizando para tener un clipping adecuado en wavsave:
  p = p/(5*rms(p));

  % Sampling frequency
  Fs = L*df;
  Fs = round(Fs);

  % Time axis
  ejet = (1/Fs)*(0:L-1);
  assert(max(ejet - linspace(0,L/Fs,L))  ...
         < (ejet(2) - ejet(1))/100, ...
         'Time axis problem.')
end

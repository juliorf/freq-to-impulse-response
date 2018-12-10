function [f, yncota] = myfft(respimp, fs)

  L = length(respimp);
 % Segun https://la.mathworks.com/help/matlab/ref/fft.html
 % sin ajustar a la potencia de dos mas cercana
  y = fft(respimp);  % Respuesta en frecuencia "original"
  ynorm = y/L; % y normalizada
  yncota = ynorm(1:(L/2)+1);   % Lado positivo del espectro
  yncota(2:end-1) = 2*yncota(2:end-1);    % x2 para mantener la energía total
  f = fs*(0:(L/2))/L; % Eje f
end

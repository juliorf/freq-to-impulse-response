function [fcota,ampcota] = myindex(f_min,f_max,f,yncota,dff)

  index_up = find(f>=f_max+dff,1,'first');
  index_low = find(f<=f_min-dff,1,'last');
  indice = find(f>=f(index_low) & f<=f(index_up));   % Indices para acotar nuestro intervalo
  fcota = f(indice);   % Eje f acotado

  ampcota = yncota(indice);
  
  assert(f_min >= f(1) & f_max <= f(end), ...
     'f_max must lie inside f')
end
function [f_i, amp_i] = interpolar_fY2(f, Y, dff, f_max, f_min)
% Interpolation of 'f' and 'Y' in a new frequency axis, starting at
% 'f_down', evenly spaced with 'f_new'. Last element of the new
% frequency axis is larger than 'f_up'.

  assert(min(f(2:end)-f(1:end-1))*1.01 ...
         > max(f(2:end)-f(1:end-1)), ...
         'Uneven frequency axis, 1/100 tolerance')

  assert(length(Y) == length(f), 'Length missmatch')

  assert(f_min >= f(1) & f_min <= f(end), ...
         'f_down must lie inside f')

%   El método de acotamiento utilizado en f2t permite que se cumpla
%   la condición del assert sin necesidad de una distancia mínima
%   (f(end) - df_new) dado que falta hacer la interpolacion y otro
%   acotamiento.
%   assert(f_max >= f(1) & f_max <= f(end) - df_new, ...
%          'f_up must lie inside [f(1), f(end) - df_new]')

    assert(f_max >= f(1) & f_max <= f(end), ...
         'f_up must lie inside [f(1), f(end)]')

  assert(f_max > f_min, ...
         'f_up must be larger than f_down')
  
  df = f(2) - f(1);
  if dff == 0 && df < 0.0248*1.01   % Parche a mano, al menos por ahora
      f_i = f;                      % El problema con este parche sería
      amp_i = Y;                    % No poder cambiar df bajo la condición
      return;
  elseif dff==0
      dff = 0.0248;
  end
  index_aux_step = dff/df;

  ii = 1;
  fdf = f_min - f(1);
  if fdf < 0
      fdf = 0;
  end
  index_aux(ii) = 1 + (fdf)/df;
  while index_aux(end) < length(f) - index_aux_step
    ii = ii + 1;
    index_aux(ii) = index_aux(ii-1) + index_aux_step;
    if f_max <= f(1)+(ii-1)*dff
      break
    end
  end
  f_i   = interp1(f, index_aux);
  amp_i = interp1(Y, index_aux);
  
% plot(f,20*log10(Y/20e-6),':o')
% hold on
% plot(f_i,20*log10(amp_i/20e-6),':*')
% hold off
% xlim([43,44])
% legend('Antes de interpolar_fY2','Despues de interpolar_fY2')
  
end

%!test
%! f = linspace(10,30,21);
%! Y = f;
%! df = f(2) - f(1);
%! df_new = df*2;
%! f_up = 20;
%! f_down = 10;
%! [f_i, a_i] = interpolar_fY(f,Y,df_new,f_up,f_down);
%! assert(f_i, [10, 12, 14, 16, 18, 20])
%! assert(f_i, a_i)

%!test
%! f = linspace(10,30,21);
%! Y = f;
%! df = f(2) - f(1);
%! df_new = df*2;
%! f_up = 20.1;
%! f_down = 10;
%! [f_i, a_i] = interpolar_fY(f,Y,df_new,f_up,f_down);
%! assert(f_i, [10, 12, 14, 16, 18, 20, 22])
%! assert(f_i, a_i)

%!test
%! f = linspace(10,30,21);
%! Y = f;
%! df = f(2) - f(1);
%! df_new = df*2;
%! f_up = 20.1;
%! f_down = 10.05;
%! [f_i, a_i] = interpolar_fY(f,Y,df_new,f_up,f_down);
%! assert(f_i, [10.05, 12.05, 14.05, 16.05, 18.05, 20.05, 22.05])
%! assert(f_i, a_i)

%!test
%! f = linspace(10,1000,1e6);
%! Y = f;
%! df = f(2) - f(1);
%! df_new = 1;
%! f_up = 400;
%! f_down = 10;
%! [f_i, a_i] = interpolar_fY(f,Y,df_new,f_up, f_down);
%! assert(f_i(end), 400, df/10)

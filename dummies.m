function [fdum,dummy] = dummies(m,f,amp,fillpad)

  if m == 0
    dummy = amp;
    fdum = f;
    return
  end
  
  fdum = f_axe_builder(m,f);
  i=1; j=1;
  if strcmp(fillpad,'zeros')
      % TODO: vectorize this loop and check performance.
      while i <= length(amp)
        dummy(j) = amp(i);  % Amplitud compleja
        for k = 1:m
          j = j+1;
          dummy(j) = 0;  % Misma amplitud
        end
        i = i+1;
        j = j+1;
      end
      
  elseif strcmp(fillpad,'repeat')
      % TODO: vectorize this loop and check performance.
      while i <= length(amp)
        dummy(j) = amp(i);  % Amplitud compleja
        for k = 1:m
          j = j+1;
          dummy(j) = amp(i);  % Misma amplitud
        end
        i = i+1;
        j = j+1;
      end
      
  elseif strcmp(fillpad,'interp')
      dummy = interp(amp,m+1);
      
  elseif strcmp(fillpad,'interp1')
      clear fdum
      fdum = f_axe_builder(m,f(1:end-1));
      fdum(end+1) = f(end);
      dummy = interp1(f,amp,fdum);
      
  elseif strcmp(fillpad,'resample')
      clear fdum
      fdum = f_axe_builder(m,f);
      dummy = resample(amp,m+1,1);
      assert(length(fdum) == length(dummy), 'Length missmatch')
      
  else
      error('Invalid parameter for fillpad')
  end
  
  i_ampstart = find(amp~=0,1,'first');
  i_dumstart = find(dummy==amp(i_ampstart),1,'first');
  dummy(1:i_dumstart-1) = 0;    % Hace cero a dummies del inicio al primer valor ~=0 de amp
  assert(length(fdum) == length(dummy),'Length missmatch')
end

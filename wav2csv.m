function [csvname,ejet,respimp] = wav2csv(filepath, nf, ...
                                            varargin)

  % TODO: para las graficas usar algo como
  % 'plots', {'plots_1', 'plots_2'}, ...
  o = parse_defaults(struct(...
    'fmin', 10, ...
    'fmax', 400, ...
    'dff', 1, ...
    'plt_w2c', false, ...
    'bpf_pass', false, ...
    'wavsave', false, ...
    'fig_visible', 'on'), ...
                     varargin{:});

  if exist('OCTAVE_VERSION', 'builtin')
    pkg load signal
  end

  assert(ischar(filepath), 'filepath: wrong type')
  assert(exist(filepath, 'file') ~= 0, [filepath, ' not found'])

  [~,filename,~] = fileparts(filepath);

  [respimp,fs] = audioread(filepath);  % Lee desde .wav f(t), fs es frec. de muestreo
  
  L = length(respimp);
  ejet = (0:L-1)/fs;  % Eje de tiempo
  assert(length(ejet) == length(respimp), 'Length missmatch')
  
  % Original impulse response plot
  if o.plt_w2c
    plots_1(ejet, respimp, ...
            filename, nf, o.fig_visible)
  end
  
  % Bandpass filter over original signal
  if o.bpf_pass
      bandpass_filter(nf,filename,ejet,respimp, ...
                      o.fig_visible,fs,o.fmin,o.fmax);
  end

  [f, yncota] = myfft(respimp, fs);
  assert(length(f) == length(yncota), 'Length missmatch')
  
  if o.plt_w2c
    plots_2(f, yncota, ...
            filename, nf, o.fig_visible, o.fmax)
  end

  % Frecuencia y amplitud acotadas
  [fcota, ampcota] = myindex(o.fmin,o.fmax,f,yncota,o.dff);
  assert(length(fcota) == length(ampcota), 'Length missmatch')
  assert(fcota(1)*1.05 >= o.fmin-o.dff | ...
         (o.fmin-o.dff <= 0 & fcota(1) == 0), ...
         'Truncated freq axis problem, low')  %*
  assert(fcota(end)*1.005 > o.fmax+o.dff, 'Truncated freq axis problem, hi')  %*

  % TODO: Check of the last element of ampcota should be purely real or
  % not.

  % Grafica resp. frec. original y acotada
  % TODO: usar
  % if ismember('plots_3', o.plots)
  if o.plt_w2c
    plots_3(f, yncota, fcota, ampcota, ...
            filename, nf, o.fig_visible)
  end

  % Cambia la diferencial de frecuencia actual por dff
  [f_i,amp_i] = interpolar_fY2(fcota,ampcota,o.dff,o.fmax,o.fmin);
  assert(abs(f_i(end) - o.fmax) ...
         < abs(mean([f_i(end),o.fmax]))/100, ...
         'Problem in interpolation')
  assert(min(f_i(2:end)-f_i(1:end-1))*1.01 ...
         > max(f_i(2:end)-f_i(1:end-1)), ...
         'Uneven frequency axis, 1/100 tolerance')

  % Compara el vector sin interpolar vs el interpolado
  if o.plt_w2c
    plots_4(fcota, ampcota, f_i, amp_i, o.dff, ...
            filename, nf, o.fig_visible)
  end

  % Crea el .csv
  csvname = strcat(nf, filesep, filename, '.csv');
  freq = f_i'; realpart = real(amp_i'); imagpart = imag(amp_i');  % Renombra las columnas
  if exist('OCTAVE_VERSION', 'builtin')
    % By now, it is not possible to use the function 'table' in
    % octave[^1], then create csv by hand
    columns = {freq', realpart', imagpart'};
    % Thanks to https://la.mathworks.com/matlabcentral/answers/263227
    fid = fopen(csvname, 'wt');
    fprintf(fid, '# Header\n');
    fclose(fid);
    % Thanks to https://stackoverflow.com/questions/14378996
    csvwrite(csvname, cell2mat(columns')', '-append');
  else
    if ~iscolumn(realpart) 
        realpart = realpart';
    end
    if ~iscolumn(imagpart) 
        imagpart = imagpart';
    end
    csvtabla = table(freq,realpart,imagpart);
    writetable(csvtabla,csvname);
  end

end

%% ---------------- Funciones locales ---------------- %%

function plots_1(ejet, respimp, ~, nf, fig_visible)

  % Grafica respuesta impulso original
  figure('visible', fig_visible)
  plot(ejet,respimp)
  title('Respuesta impulso original', ...
        'Interpreter', 'none')
  xlabel('Tiempo [s]')
  ylabel('Amplitud [Pa]')
  xlim([ejet(1),ejet(end)])
  picname = strcat(nf,filesep, ...
                   'w2c_1_resp_imp_original.png');
  print(picname, '-dpng')
end

function plots_2(f, yncota, ~, nf, fig_visible, fmax)

  figure('visible', fig_visible)
  plot(f, abs(yncota))    % Grafica el módulo de y norm. ac.
  xlim([0, fmax])
  title('y normalizada acotada', ...
        'Interpreter', 'none')
  xlabel('Frecuencia [Hz]')
  ylabel('Amplitud [Pa]')
  picname = strcat(nf, filesep, ...
                   'w2c_2.1_y_norm_acot_Pa.png');
  print(picname, '-dpng')

  figure('visible', fig_visible)
  plot(f,20*log10(abs(yncota)/20e-6)) % Grafica anterior en dB
  xlim([0, fmax])
  title('y normalizada acotada', ...
        'Interpreter', 'none')
  xlabel('Frecuencia [Hz]')
  ylabel('|Y| [dB SPL]')
  picname = strcat(nf, filesep, ...
                   'w2c_2.2_y_norm_acot_dB.png');
  print(picname, '-dpng')

end

function plots_3(f, yncota, fcota, ampcota, ~, nf, fig_visible)
  
  % Grafica resp. frec. original y acotada
  figure('visible', fig_visible)
  hold on
  plot(f,20*log10(abs(yncota)/20e-6))
  plot(fcota,20*log10(abs(ampcota)/20e-6))
  hold off
  xlim([0, 1000])
  title('y original vs y acotada', ...
        'Interpreter', 'none')
  xlabel("Frecuencia [Hz]")
  ylabel("|Y| [dB]")
  picname = strcat(nf, filesep, ...
                   'w2c_3_y_acot_sobre_original_dB.png');
  print(picname, '-dpng')
end

function plots_4(fcota, ampcota, f_i, amp_i, dff, ~, nf, fig_visible)
  % Compara el vector sin interpolar vs el interpolado
  figure('visible', fig_visible)
  plot(fcota, 20*log10(abs(ampcota)/20e-6))
  hold on
    plot(f_i, 20*log10(abs(amp_i)/20e-6))
    legend('fft directa de la senhal original', ...
         ['resampleado a ', num2str(dff), ' Hz'])
  hold off
  title('Sin interpolar vs interpolado', 'Interpreter', 'none')
  xlabel("frecuencia [Hz]")
  ylabel("|Y| [dB]")
  picname = strcat(nf, filesep, ...
                   'w2c_4.1_sin_interp1_vs_con_interp1_full.png');
  print(picname, '-dpng')
  
  figure('visible', fig_visible)
  plot(fcota, 20*log10(abs(ampcota)/20e-6),'-o')
  hold on
    plot(f_i, 20*log10(abs(amp_i)/20e-6),'--x')
    legend('fft directa de la senhal original', ...
         ['resampleado a ', num2str(dff), ' Hz'])
    xlim([106,108])
  hold off
  title('Sin interpolar vs interpolado', 'Interpreter', 'none')
  xlabel("frecuencia [Hz]")
  ylabel("|Y| [dB]")
  picname = strcat(nf, filesep, ...
                   'w2c_4.1_sin_interp1_vs_con_interp1_acotada.png');
  print(picname, '-dpng')
end

% [^1]: https://github.com/octave-de/table/issues/1
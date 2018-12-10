function p = freq2time(filepath, nf, ejet_orig, ri_orig, varargin)

  % Falta completar el "En general:"
  % En general:
  %
  % 1. Agrega ceros entre cero y el primer elemento de cada
  % vector (frecuencias, amplitudes)
  % 2. Los elementos del vector "amplitud" son módulos de los
  % números complejos del csv original
  % 3. Crea la matriz con m dummies entre cada dato original:
  % son elementos de "relleno" para alcanzar una respuesta en
  % tiempo más larga
  % 4. Concatena la matriz anterior con su complejo conjugado,
  % para [para qué?!?!]
  % 5. Aplica ifft.
  %
  % Parametros:
  %
  % m: Número de muestras entre elementos, para dummy
  %
  % plt_f2t: generar graficas para checar que los resultados sean
  % correctos
  %
  % fig_visible: mostrar las graficas en ventanas al generarlas, de
  % todos modos se generan los archivos png.

  o = parse_defaults(struct(...
    'm', 10, ...
    'fmin', 10, ...
    'fmax', 400, ...
    'plt_f2t', false, ...
    'wavsave', false, ...
    'dff', 1, ...
    'fillpad', 'interp1', ...
    'downcsv', false, ...
    'fig_visible', 'on'), ...
                     varargin{:});

  assert(ischar(filepath), 'filepath: wrong type')
  assert(exist(filepath, 'file') ~= 0, [filepath, ' not found'])
  showdff = true;   % Muestra dff en gráficas. Solo lo hacía con los wav, pero ahora con todos
  [~,filename,~] = fileparts(filepath);

  % Lee el archivo
  % En formato de Actran, el primer renglon se descarta porque tiene texto
  csv = csvread(filepath,1,0);
  
  % Vector de frecuencias
  f_csv = csv(:,1);
  aux = find(f_csv<=o.fmin,1,'last');
  aux2 = find(f_csv>=o.fmax,1,'first');
  if isempty(aux)
      aux = 1;
  end
  if isempty(aux2)
      aux2 = length(f_csv);
  end
  famp_cota = find(f_csv>=f_csv(aux) & f_csv<=f_csv(aux2));
  
  f = f_csv(famp_cota);     % Frecuencias directo de csv
  df1 = f(2)-f(1);          % Resolución
  L_freq = length(f);

  % aqui es importante construir explicitamente el vector, la notacion
  % a:b:c falla, por favor no cambiar
  ii = 1;
  padding_freq(ii) = 0;
  while padding_freq(ii) < f(1) - (1.5)*df1
    ii = ii + 1;
    padding_freq(ii) = padding_freq(ii-1) + df1;
  end
  L_pf = length(padding_freq);

  % Eje x, de cero a f final
  ejef(1:L_pf,1) = padding_freq;
  ejef(L_pf+1:L_pf+length(f),1) = f;
  
  % Construyendo el padding para arriba de f final hasta fs/2 kHz
  l_ejef = length(ejef);
  ii = 1;
  padding_freq_up(ii) = ejef(end) + df1;
  while padding_freq_up(ii) <= 22050 % fs/2
    ii = ii + 1;
    padding_freq_up(ii) = padding_freq_up(ii-1) + df1;
  end
  L_pf_up = length(padding_freq_up);
  
  % Eje x, de 0 a fs kHz
  ejef(l_ejef+1:l_ejef+L_pf_up) = padding_freq_up;

  % Crea la matriz de complejos completando entre 0 y frec. inicial
  amplitud(1:L_pf) = zeros(L_pf,1);
  amplitud(L_pf+1:L_pf+L_freq) = csv(famp_cota,2) + 1j*csv(famp_cota,3);
  amplitud(l_ejef+1:l_ejef+L_pf_up) = zeros(L_pf_up,1);
  assert(length(ejef) == length(amplitud), 'Length missmatch')
  
  % Si se va a guardar algun archivo, se crea la carpeta
  if o.plt_f2t || o.wavsave
      outdir = plots_dir(nf, filename);
  end
  
  % Baja la interpolación de la señal, si se requiere
  df1 = ejef(2)-ejef(1);
  if o.downcsv && (df1<o.dff || o.dff == 0)
      [ejef,amplitud] = interpolar_fY2(ejef, amplitud, o.dff, o.fmax, ...
                                                            o.fmin);
  elseif o.downcsv && df1>o.dff
      warning(['Frequency resolution bigger than dff, ', ...
      'this is a task for dummies.m. Interpolation aborted.'])
  end
  
  % Ventana para gráfica en frecuencia:
  f1 = 41.25;  % Limite inferior
  f2 = 43.25;    % Limite superior
  [lim_inf,lim_sup] = lim_plot(ejef,f1,f2);
  
  % Ventana para gráfica en tiempo:
  t1 = 0.3;  % Limite inferior
  t2 = 0.5;    % Limite superior
  
  % Respuesta en frecuencia original
  if o.plt_f2t
    plots_1(amplitud, ejef, filename, outdir, o.fig_visible, o.m, ...
                      o.dff, lim_inf, lim_sup, o.fmin, o.fmax, showdff)
  end
  
  % Crea la matriz con valores dummy
  [fdum,dummy] = dummies(o.m,ejef,amplitud,o.fillpad);
  assert(length(fdum) == length(dummy))

  % Después de dummies, antes de la normalización
  if o.plt_f2t
    plots_4(dummy, fdum, filename, outdir, o.fig_visible, ...
                        o.m, o.dff, f1, f2, o.fmin, o.fmax, showdff)
  end
  
  % Compara original y después de dummies como ventana
  if o.plt_f2t
    plots_5(amplitud, ejef, dummy, fdum, filename, outdir, ...
            o.fig_visible, o.m, o.dff, f1, f2, o.fmin, o.fmax)
  end
     
  dummy(2:end-1) = 0.5*dummy(2:end-1);  % Proceso inverso a la normalización
  Ndumm = length(dummy);
  % TODO: Remove this patch, from here
  if mod(length(fdum),2) ~= 1
    fdum = fdum(1:end-1);
    dummy = dummy(1:end-1);
    Ndumm = length(dummy);
  end
  % to here
  assert(mod(Ndumm,2) == 1, ...
         'By now, length of fdum must be an odd number')
     
  % Refleja la nueva matriz
  reflejada = reflex(dummy);
  reflejada(2:end) = conj(reflejada(2:end));
    
  Ldumm = 2*(Ndumm-1);
  assert(length(reflejada) == Ldumm)
  reflejada = Ldumm * reflejada;

  % Gráfica después de reflex y conjugación
  if o.plt_f2t
      plots_2(reflejada, filename, outdir, o.fig_visible, o.m, ...
                o.dff, showdff)
  end

  % Paso al dominio del tiempo
  for n = 1:length(fdum)
      df_dumm = fdum(n+1) - fdum(n);
      if df_dumm ~= 0
          break;
      end
  end
  
  [ejet,p,fs] = myifft(reflejada,df_dumm);
  assert(length(p) == Ldumm, 'Length missmatch')
  assert(length(ejet) == length(p), 'Length missmatch')
  
  % Acotando hasta el final de la respuesta original
  if (isempty(ejet_orig) && ejet(end)<=40) || ...
            (~isempty(ejet_orig) && ejet(end) <= ejet_orig(end))
      t_lim = length(ejet);
  elseif ~isempty(ejet_orig)
      t_lim = find(ejet>=ejet_orig(end),1,'first');
  elseif isempty(ejet_orig) && ejet(end)>=40
      t_lim = find(ejet<=40,1,'last');
      disp('Este caso sí ocurre')
      keyboard
  end
  
  ejet = ejet(1:t_lim);
  p = p(1:t_lim);
  assert(length(ejet) == length(p), 'Length missmatch')
  
  [lim_inf,lim_sup] = lim_plot(ejet,t1,t2);
  
  if o.plt_f2t
    plots_3(ejet, p, amplitud, ejet_orig, ri_orig, filename, outdir, ...
                o.fig_visible, o.m, o.dff, showdff, lim_inf, lim_sup)
  end

  % Respuesta impulso resultante en .wav
  if o.wavsave
      fnsave = strcat(outdir, filesep, filename, '_resultante.wav');
      audiowrite(fnsave,p,fs);
      x=[10,15,30];
      for i=1:length(x)
          fnsave = strcat(outdir, filesep, filename, '_resultante_', ...
              num2str(x(i)),'s.wav');
          tlim = find(ejet<=x(i),1,'last');
          audiowrite(fnsave,p(1:tlim),fs);
      end
  end

end

%% ---------------- Funciones locales ---------------- %%

% Magnitud y fase originales
function plots_1(amplitud, ejef, filename, nf, fig_visible, m, dff, ...
                    lim_inf, lim_sup, fmin, fmax, showdff)
  amdb = 20*log10(abs(amplitud)/20e-6);
  
  figure('visible', fig_visible)
  subplot 211
      plot(ejef,amdb)
      title(['Magnitud original de ', filename], ...
            'Interpreter', 'none')
      xlabel("Frecuencia [Hz]")
      ylabel("|Y| [dB]")
      xlim([fmin,fmax])
      legend2show(m,dff,showdff)
  subplot 212
      plot(ejef,unwrap(angle(amplitud)))
      title('Fase original','Interpreter', 'none')
      xlabel('Frecuencia [Hz]')
      ylabel('Angulo [rad]')
      xlim([fmin,fmax])
  picname = strcat(nf, filesep, ...
                   'f2t_1_resp_frec_original_full.png');
  print(picname, '-dpng')
  
  figure('visible', fig_visible)
  subplot 211
      plot(ejef,amdb,'-o','MarkerIndices',1:length(amdb))
      title(['Magnitud original de ', filename], ...
            'Interpreter', 'none')
      xlabel("Frecuencia [Hz]")
      ylabel("|Y| [dB]")
      xlim([ejef(lim_inf),ejef(lim_sup)])
      legend2show(m,dff,showdff)
  subplot 212
      plot(ejef,unwrap(angle(amplitud)),'-o','MarkerIndices',1:length(amplitud))
      title('Fase original','Interpreter', 'none')
      xlabel('Frecuencia [Hz]')
      ylabel('Angulo [rad]')
      xlim([ejef(lim_inf),ejef(lim_sup)])
  picname = strcat(nf, filesep, ...
                   'f2t_1_resp_frec_original_acotada.png');
  print(picname, '-dpng')
end

% Después de reflex
function plots_2(reflejada, ~, nf, fig_visible, m, dff, showdff)
                    
  figure('visible', fig_visible)
  plot(abs(reflejada))
  xlabel("numero de entrada")
  ylabel("|Y| [Pa]")
  title('Despues de reflex, modulo', 'Interpreter', 'none')
  legend2show(m,dff,showdff)
  picname = strcat(nf, filesep, 'f2t_2.1_reflejada_abs.png');
  print(picname, '-dpng')

  figure('visible', fig_visible)
  plot(real(reflejada))
  xlabel("numero de entrada")
  ylabel("real(Y) [Pa]")
  title('Despues de reflex, parte real', 'Interpreter', 'none')
  legend2show(m,dff,showdff)
  picname = strcat(nf, filesep, 'f2t_2.2_reflejada_real.png');
  print(picname, '-dpng')

  figure('visible', fig_visible)
  plot(imag(reflejada))
  xlabel("numero de entrada")
  ylabel("imag(Y) [Pa]")
  title('Despues de reflex, parte imaginaria','Interpreter', 'none')
  legend2show(m,dff,showdff)
  picname = strcat(nf, filesep, 'f2t_2.3_reflejada_imag.png');
  print(picname, '-dpng')
end

% Respuesta impulso original y resultante
function plots_3(ejet, p, ~, ejet_orig, ri_orig, filename, ...
                    nf, fig_visible, m, dff, showdff, lim_inf, lim_sup)
    figure('visible', fig_visible)
    plot(ejet,p)
    title(['R.I. resultante de ', filename], 'Interpreter', 'none')
    xlabel('Tiempo [s]')
    ylabel('Amplitud [Pa]')
    legend2show(m,dff,showdff)
    legend('Location','southeast','Orientation','horizontal')
    picname = strcat(nf, filesep, ...
                   'f2t_3.1_Respuesta_impulso_full.png');
    print(picname, '-dpng')

    figure('visible', fig_visible)
    plot(ejet,p,'-o')
    title(['R.I. resultante de ', filename], 'Interpreter', 'none')
    xlabel('Tiempo [s]')
    ylabel('Amplitud [Pa]')
    legend2show(m,dff,showdff)
    legend('Location','southeast','Orientation','horizontal')
    xlim([ejet(lim_inf),ejet(lim_sup)])
    picname = strcat(nf, filesep, ...
                   'f2t_3.1_Respuesta_impulso_acotada.png');
    print(picname, '-dpng')

    % Compara wav vs resultante
    figure('visible', fig_visible)
    rimax = ri_orig/max(ri_orig);
    plot(ejet_orig,rimax/max(abs(rimax)),ejet,p/max(abs(p)))
    title(['Respuesta impulso de ',filename],'Interpreter', 'none')
    xlabel('Tiempo [s]')
    ylabel('Amplitud [Pa]')
    legend('Original', 'Resultante')
    legend('Location','southeast','Orientation','horizontal')
    picname = strcat(nf, filesep, ...
                     'f2t_3.2_RI_original_vs_result_full.png');
    print(picname, '-dpng')

    figure('visible', fig_visible)
    rimax = ri_orig/max(ri_orig);
    plot(ejet_orig,rimax/max(abs(rimax)),'-o',ejet,p/max(abs(p)),':x')
    title(['Respuesta impulso de ',filename],'Interpreter', 'none')
    xlabel('Tiempo [s]')
    ylabel('Amplitud [Pa]')
    legend('Original', 'Resultante')
    legend('Location','southeast','Orientation','horizontal')
    xlim([ejet(lim_inf),ejet(lim_sup)])
    picname = strcat(nf, filesep, ...
                     'f2t_3.2_RI_original_vs_result_acotada.png');
    print(picname, '-dpng')
end

% Respuesta en frecuencia con dummies
function plots_4(dummy, fdum, filename, nf, fig_visible, m, dff, ...
                    f1, f2, fmin, fmax, showdff)
  dumdb = 20*log10(abs(dummy)/20e-6);
  [lim_inf,lim_sup] = lim_plot(fdum,f1,f2);
  
  figure('visible', fig_visible)
  subplot 211
      plot(fdum,dumdb)
      title(['Magnitud resultante de dummies para ', filename], ...
            'Interpreter', 'none')
      xlabel("Frecuencia [Hz]")
      ylabel("|Y| [dB]")
      xlim([fmin,fmax])
      legend2show(m,dff,showdff)
  subplot 212
      plot(fdum,unwrap(angle(dummy)))
      title('Fase resultante','Interpreter', 'none')
      xlabel('Frecuencia [Hz]')
      ylabel('Angulo [rad]')
      xlim([fmin,fmax])
  picname = strcat(nf, filesep, 'f2t_4_resp_frec_dummies_full.png');
  print(picname, '-dpng')
  
  figure('visible', fig_visible)
  subplot 211
      plot(fdum,dumdb,'-o','MarkerIndices',1:length(dumdb))
      title(['Magnitud resultante de dummies para ', ...
             filename], ...
            'Interpreter', 'none')
      xlabel("Frecuencia [Hz]")
      ylabel("|Y| [dB]")
      xlim([fdum(lim_inf),fdum(lim_sup)]) % Intervalo definido por f1 y f2
      legend2show(m,dff,showdff)
  subplot 212
      plot(fdum,unwrap(angle(dummy)),'-o','MarkerIndices',1:length(dumdb))
      title('Fase resultante','Interpreter', 'none')
      xlabel('Frecuencia [Hz]')
      ylabel('Angulo [rad]')
      xlim([fdum(lim_inf),fdum(lim_sup)]) % Intervalo definido por f1 y f2
  picname = strcat(nf, filesep, 'f2t_4_resp_frec_dummies_acotada.png');
  print(picname, '-dpng')
end

% Respuesta en frecuencia original vs. dummies
function plots_5(amplitud, ejef, dummy, fdum, filename, outdir, ...
                 fig_visible, ~, ~, f1, f2, fmin, fmax)
  amdb = 20*log10(abs(amplitud)/20e-6);
  dumdb = 20*log10(abs(dummy)/20e-6);
  [lim_inf,lim_sup] = lim_plot(fdum,f1,f2);
  
  figure('visible', fig_visible)
  subplot 211
      plot(ejef,amdb,'LineWidth',2)
      hold on
        plot(fdum,dumdb)
      hold off
      title(['Magnitud original vs. dummies para ', filename], 'Interpreter', 'none')
      xlabel('Frecuencia [Hz]')
      ylabel('|Y| [dB]')
      xlim([fmin,fmax])
      legend('Original','Dummies')
  subplot 212
      plot(ejef,unwrap(angle(amplitud)),'LineWidth',2)
      hold on
        plot(fdum,unwrap(angle(dummy)))
      hold off
      title('Fase')
      xlabel('Frecuencia [Hz]')
      ylabel('Angulo [rad]')
      xlim([fmin,fmax])
      legend('Original','Dummies')
  picname = strcat(outdir, filesep, ...
                   'f2t_5_resp_frec_original_vs_dummies_full.png');
  print(picname, '-dpng')
  
  figure('visible', fig_visible)
  subplot 211
      plot(ejef,amdb,'-o','MarkerIndices',1:length(amdb),'LineWidth',2)
      hold on
        plot(fdum,dumdb,'--x','MarkerIndices',1:length(fdum))
      hold off
      title(['Magnitud original vs. dummies para ', filename], 'Interpreter', 'none')
      xlabel('Frecuencia [Hz]')
      ylabel('|Y| [dB]')
      xlim([fdum(lim_inf),fdum(lim_sup)]) % Intervalo definido por f1 y f2
      legend('Original','Dummies')
  subplot 212
      plot(ejef,unwrap(angle(amplitud)),'-o','MarkerIndices', ...
          1:length(amplitud),'LineWidth',2)
      hold on
        plot(fdum,unwrap(angle(dummy)),'--x','MarkerIndices', ...
            1:length(fdum))
      hold off
      title('Fase')
      xlabel('Frecuencia [Hz]')
      ylabel('Angulo [rad]')
      xlim([fdum(lim_inf),fdum(lim_sup)]) % Intervalo definido por f1 y f2
      legend('Original','Dummies')
  picname = strcat(outdir, filesep, ...
                   'f2t_5_resp_frec_original_vs_dummies_acotada.png');
  print(picname, '-dpng')
end

function [lim_inf,lim_sup] = lim_plot(ejex,x1,x2)
    lim_inf = find(ejex<=x1,1,'last');
    lim_sup = find(ejex>=x2,1,'first');
    if (x1<ejex(1) || x1>ejex(end)) && (x2<ejex(1) || x2>ejex(end))
        lim_inf = 1;
        lim_sup = length(ejex);
    end
end

function legend2show(m,dff,showdff)
    if showdff
        legend(['m=',num2str(m),', dff=',num2str(dff)])
    else
        legend(['m=',num2str(m)])
    end
end
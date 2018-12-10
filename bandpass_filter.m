function bandpass_filter(nf,filename,t,amp,fig_visible,fs,fmin,fmax)

    % Filter design
    d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2,C',...
                    117,9.99,10,400,400.01,10000);
    d.Stopband1Constrained = true;
    d.Astop1 = 60;
    d.Stopband2Constrained = true;
    d.Astop2 = 60;

    Hd = design(d,'equiripple');    % Passband type
    h = fvtool(Hd);    % Filter design plot
    legend(h,'Orden del filtro: 117')
    picname = strcat(nf, filesep, ...
                     'w2c_BPF_filter.png');
    print(picname, '-dpng')
%     h.visible = 'off';

    % Filter apply
    s_fil = filter(Hd,amp);
    
    % Creates .wav from results
%     if wavsave     % TODO: ajustar para pedir confirmacion
      fnsave = strcat(nf,'\',filename,'_filtrado_resultante.wav');
      audiowrite(fnsave,s_fil,fs);
%     end

    % Results plot, time domain
    figure('visible',fig_visible)
    
    y_max = max(abs(amp));
    subplot(2,1,1)
        plot(t,amp)
        title(['R.I. original de ', filename], ...
              'Interpreter', 'none')
        xlabel('Tiempo [s]')
        ylabel('Amplitud [Pa]')
        xlim([0,t(end)])
        ylim([-y_max, y_max])

    limx = find(t==10);
    y_max = max(abs(s_fil));
    subplot(2,1,2)
        plot(t,s_fil)
        title('R.I. filtrada','Interpreter','none')
        xlabel('Tiempo [s]')
        ylabel('Amplitud [Pa]')
        xlim([0,t(limx)])
        ylim([-y_max, y_max])

    picname = strcat(nf, filesep, ...
                     'w2c_BPF_1_R.I._original_vs_filtrada.png');
    print(picname, '-dpng')
    
    % fft to amp (original)
    ampfft = fft(amp);  % original signal's fft
    l_ampfft = length(ampfft);
    ampfnorm = abs(ampfft/l_ampfft);
    ampfncut = ampfnorm(1:l_ampfft/2+1);
    ampfncut(2:end-1) = 2*ampfncut(2:end-1);
    
    sf_fft = fft(s_fil);   % filtered signal's fft
    sf_fnorm = abs(sf_fft/l_ampfft);
    sf_fncut = sf_fnorm(1:l_ampfft/2+1);
    sf_fncut(2:end-1) = 2*sf_fncut(2:end-1);
    
    f = fs*(0:(l_ampfft/2))/l_ampfft;   % Freq axis
    f_smaller = find(f<fmin);  % Working interval
    f_higher = find(f>fmax);
    
    figure('visible',fig_visible)
    subplot(2,1,1)
        plot(f,ampfncut)
        title(['Resp. en frec. original de ', filename], ...
              'Interpreter', 'none')
        xlabel('Frecuencia [Hz]')
        xlim([f(f_smaller(1)),f(f_higher(1))+500])
        ylabel('Amplitud [Pa]')
    subplot(2,1,2)
        plot(f,sf_fncut)
        title('Espectro de frecuencia después del filtro','Interpreter','none')
        xlabel('Frecuencia [Hz]')
        xlim([f(f_smaller(1)),f(f_higher(1))+500])
        ylabel('Amplitud [Pa]')

    picname = strcat(nf, filesep, ...
                     'w2c_BPF_2_R.I._frecuenciaoriginales.png');
    print(picname, '-dpng')
    
    % "Cutting off" the values out of our range
%     ampfncut(f_smaller) = 0;
%     ampfncut(f_higher) = 0;
%     sf_fncut(f_smaller) = 0;
%     sf_fncut(f_higher) = 0;
%     
%     % Going to dB SPL
%     amp_db = 20*log10(abs(ampfncut)/20e-6);
%     sf_db = 20*log10(abs(sf_fncut)/20e-6);
%     
%     figure('visible',fig_visible)
%     subplot(2,1,1)
%         plot(f,amp_db)
%         title(['Resp. en frec. original de ', filename], ...
%               'Interpreter', 'none')
%         xlabel('Frecuencia [Hz]')
%         xlim([0,f(f_higher(1))+500])
%         ylabel('Amplitud [Pa]')
%     subplot(2,1,2)
%         plot(f,sf_db)
%         title('Espectro de frecuencia después del filtro','Interpreter','none')
%         xlabel('Frecuencia [Hz]')
%         xlim([0,f(f_higher(1))+500])
%         ylabel('Amplitud [Pa]')
% 
%     picname = strcat(nf, filesep, ...
%                      'w2c_BPF_2_R.I._frecuencia.png');
%     print(picname, '-dpng')

end
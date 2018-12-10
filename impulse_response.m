clc; close all
clear

% dummies = [40,4,0];  % Numero de dummies/Factor de incremento en interp
new_df = [0.1,1,0];    % Nueva df. df original: 0.0248
fmin = 10;       % [Hz] 
fmax = 400;      % [Hz]
plot_w2c = true;
plot_f2t = true;
bpf_pass = false;   % If true, filtrates and plots results.
                    % Otherwise, skips bandpass_filter function
wavsave = true;     % Creates a .wav file from resultant IR
downcsv = true;     % If true, turns csv's df into dff

% Method to add dummies:
    % 'repeat' for dummies to repeat amplitude
    % 'zeros' for amplitude equal to zero
    % 'interp' for interp function
    % 'interp1' for interp1 function
    % 'resample' for resample function
methods = ["repeat", "zeros", "interp", "interp1", "resample"];  

nf = get_outdir();  % Crea carpeta de trabajo
if ~exist(nf,'dir')
    mkdir(nf);
end

warningsbye()
for j=1:length(new_df)  
    dff = new_df(j);
    if dff == 0.1
        m=4;
    elseif dff == 1
        m=40;
    elseif dff == 0
        m=0;
    end
    
    disp(['Current: m = ', num2str(m), ', dff = ', num2str(dff)])
    dir_name = strcat('m', string(m), '_dff', string(dff));
    curdir = plots_dir(nf,char(dir_name));

    wav_filepath = 'Resp_Imp_Der_Abajo.wav';
    % Pasa el wav a csv
    [csv_filepath, ejet_orig, ri_orig] = ...
        wav2csv(wav_filepath, curdir, ...
                 'plt_w2c', plot_w2c, ...
                 'bpf_pass', bpf_pass, ...
                 'wavsave', wavsave, ...
                 'fig_visible', 'off', ...
                 'dff', dff, ...
                 'fmin', fmin, ...
                 'fmax', fmax);

    for k=1:length(methods)
        fillpad = methods(k);
        mthdir = plots_dir(curdir,char(fillpad));
        disp(['- Method: ', char(fillpad)])

        if dff == 0
            csv_filepath = strcat(curdir,filesep,'Resp_Imp_Der_Abajo.csv');
            csv_filepath = char(csv_filepath);
            disp('-- Now: Resp_Imp_Der_Abajo.csv')
            freq2time(csv_filepath, mthdir, ejet_orig, ri_orig, ...
                      'm', m, ...
                      'fmin', fmin, ...
                      'fmax', fmax, ...
                      'plt_f2t', plot_f2t, ...
                      'wavsave', wavsave, ...
                      'dff', dff, ...
                      'fillpad', fillpad, ...
                      'fig_visible', 'off');
        end

        % csv de Actran
        if dff==1   % Executes exclusively to get df aprox. to 0.0248
            csv_filepath = 'actran_daniel_micro2.csv';
            disp(['-- Now: ', char(csv_filepath)])

            freq2time(csv_filepath, mthdir, ejet_orig, ri_orig, ...
                      'm', m, ...
                      'fmin', fmin, ...
                      'fmax', fmax, ...
                      'plt_f2t', plot_f2t, ...
                      'wavsave', wavsave, ...
                      'dff', dff, ...
                      'fillpad', fillpad, ...
                      'downcsv', false, ...
                      'fig_visible', 'off');
        end

        % csv de FEniCS
        csv_filepath = 'micro2Daniel_0_400xx.csv';
        disp(['-- Now: ', char(csv_filepath)])

        freq2time(csv_filepath, mthdir, ejet_orig, ri_orig, ...
                  'm', m, ...
                  'fmin', fmin, ...
                  'fmax', fmax, ...
                  'plt_f2t', plot_f2t, ...
                  'wavsave', wavsave, ...
                  'dff', dff, ...
                  'fillpad', fillpad, ...
                  'downcsv', downcsv, ...
                  'fig_visible', 'off');

        if m==0
            break; % Runs only once to not repeat results
        end
    end
end

function warningsbye
    warning('off','MATLAB:colon:nonIntegerIndex')
    warning('off','MATLAB:audiovideo:audiowrite:dataClipped')
end
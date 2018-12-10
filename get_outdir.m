function outdir = get_outdir()

  % This is not optimal, see https://stackoverflow.com/questions/45461173
  % TODO: remove this and use &&, even in matlab. Remove code below [^1]
  if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:possible-matlab-short-circuit-operator');
  end

  [~,b]=system('git rev-parse HEAD');

  short_commit_hash = b(1:8);
  dt = datestr(now,'yyyy-mm-dd_HH.MM.SS');
  outdir_commit = strcat('output_commit_', dt, '_', short_commit_hash);

  % Default value
  outdir = outdir_commit;

  files = dir(pwd);
  for i=1:length(files)
    if exist(files(i).name, 'dir')==7 ...
          & regexp(files(i).name, short_commit_hash)
      % override default if already there
      outdir = 'output_dirty';
      break
    end
  end

  disp(['Writing results in ''', outdir, ''''])

  % [^1]: TODO: remove this.
  if exist('OCTAVE_VERSION', 'builtin')
    warning('on', 'Octave:possible-matlab-short-circuit-operator');
  end

end

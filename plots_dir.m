function outdir = plots_dir(path, dirname)
  outdir = [path, filesep, dirname];
  if ~exist(outdir,'dir')
    mkdir(outdir);
  end
end

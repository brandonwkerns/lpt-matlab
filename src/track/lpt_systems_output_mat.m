function lpt_systems_output_mat(systems, matfile)

    eval(['!rm -f ',matfile])

    %CE = objects;
    TIMECLUSTERS = systems;

    %% .mat file output
    %% If "TIMECLUSTERS" is > 2 GB, need to break it up.

    sum_of_mb = 0.0;
    for tt = 1:numel(TIMECLUSTERS)
      TCtemp = TIMECLUSTERS(tt);
      stats = whos('TCtemp');
      sum_of_mb = sum_of_mb + stats.bytes/1024/1024;
    end

    disp(['Total TIMECLUSTERS size: ',num2str(sum_of_mb),' MB.'])

    n_breakup = ceil(sum_of_mb / 2000.0);

    disp(['Will break it up in to around ', n_breakup, ' structs like TIMECLUSTERS, TIMECLUSTERS2, TIMECLUSTERS3, ect.'])

    nn = 1;
    more_to_break_up = 1;

    while more_to_break_up
    %for nn = 1:n_breakup

      break_up_point = -1;
      sum_of_mb = 0.0;

      if nn == 1
        tc0 = 'TIMECLUSTERS';
      else
        tc0 = ['TIMECLUSTERS', num2str(nn)];
      end

      tc1 = ['TIMECLUSTERS', num2str(nn+1)];

      eval(['TC0 = ', tc0, ';'])

      for tt = 1:numel(TC0)
      	TCtemp = TC0(tt);
      	stats = whos('TCtemp');
      	sum_of_mb = sum_of_mb + stats.bytes/1024/1024;
        if (sum_of_mb > 2000.0)
          break_up_point = tt-1;
          break
        end
      end

      if break_up_point < 0
        eval(['fout.',tc0,' = TC0']) ;
        more_to_break_up = 0;
      else
      	disp(['Data larger than 2 GB! Broken up in to ',tc0, ' and ',tc1,'.'])
      	eval(['fout.',tc0,' = TC0(1:break_up_point);'])
      	eval(['fout.',tc1,' = TC0(break_up_point+1:end);'])
      	eval([tc1, ' = TC0(break_up_point+1:end);'])
        nn = nn + 1;
      end
    end

    %fout.grid = CE.grid ;
    disp(matfile)
    eval(['save ', matfile, ' -struct fout'])

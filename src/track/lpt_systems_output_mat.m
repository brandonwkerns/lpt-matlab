function lpt_systems_output_mat(systems, objects, matfile)

    eval(['!rm -f ',matfile])

    CE = objects;
    TIMECLUSTERS = systems;
    
		   % .mat file output
		   % If "TIMECLUSTERS" is > 2 GB, need to break it up.

    sum_of_mb = 0.0;
    for tt = 1:numel(TIMECLUSTERS)
      TCtemp = TIMECLUSTERS(tt);
      stats = whos('TCtemp');
      sum_of_mb = sum_of_mb + stats.bytes/1024/1024;
    end
    
    disp(['Total TIMECLUSTERS size: ',num2str(sum_of_mb),' MB.'])

    n_breakup = ceil(sum_of_mb / 2000.0);
    
    disp(['Will break it up in to ',n_breakup, ' structs like TIMECLUSTERS, TIMECLUSTERS2, TIMECLUSTERS3, ect.'])

    for nn = 1:n_breakup
    
      break_up_point = -1;
      sum_of_mb = 0.0;

      if nn == 1
	tc0 = 'TIMECLUSTERS';
      else
	tc0 = ['TIMECLUSTERS', num2str(nn)];
      end

      tc1 = ['TIMECLUSTERS', num2str(nn+1)];

      eval(['TC0 = ', tc0, ';']) 
      %eval(['TC1 = ', tc1, ';']) 
      
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
      else
	disp(['Data larger than 2 GB! Broken up in to ',tc0, ' and ',tc1,'.'])
	eval(['fout.',tc0,' = TC0(1:break_up_point);'])
	eval(['fout.',tc1,' = TC0(break_up_point+1:end);'])
	eval([tc1, ' = TC0(break_up_point+1:end);']) 	
      end
    
    end

    %{
    if isfield(fout, 'TIMECLUSTERS2')
      
      TIMECLUSTERS2 = fout.TIMECLUSTERS2;
      
      break_up_point = -1;
      sum_of_mb = 0.0;
      for tt = 1:numel(TIMECLUSTERS2)
	TCtemp = TIMECLUSTERS2(tt);
	stats = whos('TCtemp');
	sum_of_mb = sum_of_mb + stats.bytes/1024/1024;
	if (sum_of_mb > 2000.0)
	  break_up_point = tt-1;
	  break
	end
      end
      
      if break_up_point < 0
	fout.TIMECLUSTERS2 = TIMECLUSTERS2 ;
      else
	disp(['Data larger than 2 GB! Broken up in to TIMECLUSTERS2 and TIMECLUSTERS3.'])
	fout.TIMECLUSTERS2 = TIMECLUSTERS2(1:break_up_point) ;
	fout.TIMECLUSTERS3 = TIMECLUSTERS2(break_up_point+1:end) ;
      end
      
    end
    
    
    if isfield(fout, 'TIMECLUSTERS3')
      
      TIMECLUSTERS3 = fout.TIMECLUSTERS3;
      
      break_up_point = -1;
      sum_of_mb = 0.0;
      for tt = 1:numel(TIMECLUSTERS3)
	TCtemp = TIMECLUSTERS3(tt);
	stats = whos('TCtemp');
	sum_of_mb = sum_of_mb + stats.bytes/1024/1024;
	if (sum_of_mb > 2000.0)
	  break_up_point = tt-1;
	  break
	end
      end
      
      if break_up_point < 0
	fout.TIMECLUSTERS3 = TIMECLUSTERS3 ;
      else
	disp('Data larger than 2 GB! Broken up in to TIMECLUSTERS3 and TIMECLUSTERS4.')
	fout.TIMECLUSTERS3 = TIMECLUSTERS3(1:break_up_point) ;
	fout.TIMECLUSTERS4 = TIMECLUSTERS3(break_up_point+1:end) ;
      end
      
    end
    
    %}
    
    fout.grid = CE.grid ;
    disp(matfile)
    eval(['save ', matfile, ' -struct fout'])
    

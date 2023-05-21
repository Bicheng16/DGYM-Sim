function Run_Finished(CPU_time_start, Chop_time_counter, flag_run_finish_ok)

    %% CPU time print
    CPU_time_end = clock;
    CPU_time = etime(CPU_time_end, CPU_time_start);
    disp(  strcat('CPU time: ', num2str(CPU_time), ' sec')    )
    disp(  strcat('Time chops performed: ', num2str(Chop_time_counter))    )

    %% Final run messages
    if flag_run_finish_ok==false
        disp('Run Finished with ERRORS')
    else
        disp('Run Finished OK')
    end

end
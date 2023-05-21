function writeReport(output, rsm_fd, rst, fluid, init, current_sim_t)

    % This function is to write a summary report for initial and each
    % timestep 
       
    file_path= strcat(pwd, '\', rsm_fd);    
    
    % Creat rpt file
    if current_sim_t == 0
        delete(fullfile(file_path, '\*.rpt'));        
        o_fid = fopen(fullfile(file_path, output{1}), 'w');
    else
        o_fid = fopen(fullfile(file_path, output{1}), 'a+');        
    end
          
    % Write title for the simulation report 
    fprintf(o_fid, '%88s\n', '****************************************************************************************');
    
    if current_sim_t == 0
        fprintf(o_fid, '%23s\n\n', 'INITIAL FLUIDS IN PLACE');
        fprintf(o_fid, '%37s\n\n', 'TOTAL FLUIDS IN PLACE FOR WHOLE FIELD');

        % Moles of each component in place
        for iComp = 1 : fluid.nComps+1
            % Get composition name 
            if iComp == fluid.nComps+1
                comp_name = 'WATER';
            else
                comp_name = fluid.comp_name{iComp};
            end

            fprintf(o_fid, '%16s %6s %4s %13.7e %11s\n', ...
            '  TOTAL INITIAL ', comp_name, ' =  ', rst.COMPIP(iComp), '    lb-mole' );                 
        end

        % Fluid in place in each phase 
        free_gas_frac = rst.FGIP / (rst.FGIP + rst.FAGIP) * 100; 
        fprintf(o_fid, '\n%36s %13.7e %7s\n', '  TOTAL INITIAL AQUEOUS  PHASE   =  ', rst.FWIP, '    STB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL INITIAL  LIQUID  PHASE   =  ', rst.FOIP, '    STB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL INITIAL   VAPOR  PHASE   =  ', rst.FGIP, '    MCF');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL INITIAL ADSORBED  PHASE  =  ', rst.FAGIP, '    MCF');
        fprintf(o_fid, '%36s %5f %2s \n\n', '  FREE GAS FRACTION  =  ', free_gas_frac, '   %' );
        
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL INITIAL PORE VOLUME      =  ', rst.PV,   '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL INITIAL HC PORE VOLUME   =  ', rst.HCPV, '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL INITIAL WATER            =  ', rst.VWat, '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL INITIAL LIQUID OIL       =  ', rst.VOil, '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n\n', '  TOTAL INITIAL FREE GAS VOLUME  =  ', rst.VGas, '    RVB');        
        
    else
        fprintf(o_fid, '\n%7s %7f %6s\n\n', 'TIME = ', current_sim_t, '  DAYS');
        fprintf(o_fid, '%37s\n\n', 'TOTAL FLUIDS IN PLACE FOR WHOLE FIELD');

        % Moles of each component in place
        for iComp = 1 : fluid.nComps+1
            % Get composition name 
            if iComp == fluid.nComps+1
                comp_name = 'WATER';
            else
                comp_name = fluid.comp_name{iComp};
            end

            fprintf(o_fid, '%16s %6s %4s %13.7e %11s\n', ...
            '  TOTAL PRESENT ', comp_name, ' =  ', rst.COMPIP(iComp), '    lb-mole' );                 
        end

        % Fluid in place in each phase 
        free_gas_frac = rst.FGIP / (rst.FGIP + rst.FAGIP) * 100;
        fprintf(o_fid, '\n%36s %13.7e %7s\n', '  TOTAL PRESENT AQUEOUS   PHASE  =  ', rst.FWIP, '    STB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL PRESENT  LIQUID   PHASE  =  ', rst.FOIP, '    STB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL PRESENT   VAPOR   PHASE  =  ', rst.FGIP, '    MCF');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL PRESENT ADSORBED  PHASE  =  ', rst.FAGIP, '    MCF');
        fprintf(o_fid, '%36s %5f %7s \n\n', '  FREE GAS FRACTION  =  ', free_gas_frac, '   %' );
        
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL PRESENT PORE VOLUME      =  ', rst.PV,   '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL PRESENT HC PORE VOLUME   =  ', rst.HCPV, '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL PRESENT WATER            =  ', rst.VWat, '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  TOTAL PRESENT LIQUID OIL       =  ', rst.VOil, '    RVB');
        fprintf(o_fid, '%36s %13.7e %7s\n\n', '  TOTAL PRESENT FREE GAS VOLUME  =  ', rst.VGas, '    RVB');
        
        % Cumulative fluid production
        fprintf(o_fid, '%44s\n\n', 'CUMULATIVE PRODUCTION OF PHASE AND COMPONENT');
        
        fprintf(o_fid, '%36s %13.7e %7s\n', '  CUMULATIVE WATER PRODUCTION    =  ', rst.FWPT, '    STB');
        fprintf(o_fid, '%36s %13.7e %7s\n', '  CUMULATIVE LIQUID PRODUCTION   =  ', rst.FOPT, '    STB');
        fprintf(o_fid, '%36s %13.7e %7s\n\n', '  CUMULATIVE GAS PRODUCTION      =  ', rst.FGPT, '    MCF');
        
        % Calculate recovery of each phase 
        recovery_wat = rst.FWPT / init.OWIP * 100;               % Water recovery  
        recovery_oil = rst.FOPT / init.OOIP * 100;               % Oil recovery
        recovery_gas = rst.FGPT / init.OGIP * 100;               % Gas recovery 
        recovery_wat(isnan(recovery_wat) | isinf(recovery_wat)) = 0; 
        recovery_oil(isnan(recovery_oil) | isinf(recovery_oil)) = 0; 
        recovery_gas(isnan(recovery_gas) | isinf(recovery_gas)) = 0;        

        fprintf(o_fid, '%42s\n\n', 'CUMULATIVE RECOVERY OF PHASE AND COMPONENT');
        fprintf(o_fid, '%36s %5f %2s \n', '  CUMULATIVE WATER RECOVERY      =  ', recovery_wat, ' %' );
        fprintf(o_fid, '%36s %5f %2s \n', '  CUMULATIVE OIL RECOVERY        =  ', recovery_oil, ' %' );
        fprintf(o_fid, '%36s %5f %2s \n', '  CUMULATIVE GAS RECOVERY        =  ', recovery_gas, ' %' );
                
    end
    
    fprintf(o_fid, '\n');
                           
    fclose(o_fid);
    
    % 
    % Creat rsm file
    %
    if current_sim_t == 0
        delete(fullfile(file_path, '\*.rsm'));
        o_fid = fopen(fullfile(file_path, output{2}), 'w');
        iname = strsplit(output{2}, '.');
        iiname = char(iname(1));
        fprintf(o_fid, '%16s\n %15s\n', ' SUMMARY OF RUN ', iiname);
        fprintf(o_fid, '%175s\n', '------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
        fprintf(o_fid, '%175s\n', 'TIME              TIME STEP     Newtons   FPR                 FOSAT            FGSAT            FWSAT            FOPR             FWPR             FGPR             PROD BHP  ');
        fprintf(o_fid, '%175s\n', 'DAYS              DAYS                    PSIA                                                                   STB/DAY          STB/DAY          MSCF/DAY         PSIA      ');
        fprintf(o_fid, '%175s\n', '------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
        fclose(o_fid);         
    else
        o_fid = fopen(fullfile(file_path, output{2}), 'a+');        

        % Production table
        space = '     ';
        years = current_sim_t / 365;
        fprintf(o_fid, '%-12.3f %5s %-8.3f %5s %2u %5s %-12.7e %5s %-9.4e %5s %-9.4e %5s %-9.4e %5s %-9.4e %5s %-9.4e %5s %-9.4e %5s %-9.4e', ...
        current_sim_t, space, rst.dt, space, rst.newtons, space, rst.P_res_HCPVweigthed, space, rst.So_PVweigthed, ...
        space, rst.Sg_PVweigthed, space, rst.Sw_PVweigthed, space, rst.FOPR, space, rst.FWPR, ...
        space, rst.FGPR, space, rst.well_P(1) );

        fprintf(o_fid, '\n');                           
        fclose(o_fid);    

    end
    
    % Assume the second well is an injector
    nWells = size(rst.well_P, 1);
    if nWells > 1
        if current_sim_t == 0
            o_fid = fopen(fullfile(file_path, output{3}), 'w');
            iname = strsplit(output{3}, '.');
            iiname = char(iname(1));
            fprintf(o_fid, '%16s\n %15s\n', ' SUMMARY OF RUN ', iiname);
            fprintf(o_fid, '%82s\n', '---------------------------------------------------------------------------------');
            fprintf(o_fid, '%82s\n', 'TIME               FOIR             FWIR            FGIR              INJE BHP   ');
            fprintf(o_fid, '%82s\n', 'DAYS               STB/DAY          STB/DAY         MSCF/DAY          PSIA       ');
            fprintf(o_fid, '%82s\n', '---------------------------------------------------------------------------------');
            fclose(o_fid);         
        else
            o_fid = fopen(fullfile(file_path, output{3}), 'a+');        

            % Injection table 
            fprintf(o_fid, '%-12.3f %5s %-9.4e %5s %-9.4e %5s %-9.4e %5s %-9.4e', ...
            current_sim_t, space, rst.FOIR, space, rst.FWIR, space, rst.FGIR, space, rst.well_P(2) );

            fprintf(o_fid, '\n');                           
            fclose(o_fid);         
        end        
                
    end
    
    
end
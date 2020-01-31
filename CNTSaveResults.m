%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        CNTSaveResults.m 
% writes results for the carbon ananotube transistor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results_File = ['mos_results-',date];
% save(Results_File, 'XI', 'E_sub', 'n', 'Ne_sub', 'Vd_bias', 'Vg_bias', 'phi_g', 'Ie', 'Ng_step', 'Nd_step', 'Lg', 't_ins', 'Egh1',...
%                     'En_old', 'band_old', 'Gn_old', 'Gp_old', 'myEE', 'myTr', 'Trraspr');
% equivalent to what is written in 'mos_results'

fprintf(Err_File_H,'\n\n*********************** Simulation Complete ***********************\n');
fprintf(Err_File_H,'\nSimulation start time: %d/%d/%d - %d:%d:%d \n', ...
    Start_time(1),Start_time(2),Start_time(3),Start_time(4),Start_time(5),Start_time(6));

fprintf(1,'\n\n*********************** Simulation Complete ***********************\n');
fprintf(1,'\nSimulation start time: %d/%d/%d - %d:%d:%d \n', ...
    Start_time(1),Start_time(2),Start_time(3),Start_time(4),Start_time(5),Start_time(6));
End_time = fix(clock);

fprintf(Err_File_H,'\nSimulation end time: %d/%d/%d - %d:%d:%d \n', ...
    End_time(1),End_time(2),End_time(3),End_time(4),End_time(5),End_time(6));
fprintf(1,'\nSimulation end time: %d/%d/%d - %d:%d:%d \n', ...
    End_time(1),End_time(2),End_time(3),End_time(4),End_time(5),End_time(6));

fclose(Err_File_H);
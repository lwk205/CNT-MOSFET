%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        CNTInitialWrite.m 
% opens file and inital output for the carbon ananotube transistor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

File_name = ['Error_File-',date,'.txt'];        % Open error file
Err_File_H = fopen(File_name,'w');
Start_time = fix(clock);

fprintf(1,'\n\n Simulation start time: %d/%d/%d - %d:%d:%d \n', ...
    Start_time(1),Start_time(2),Start_time(3),Start_time(4),Start_time(5),Start_time(6));
fprintf(Err_File_H,'\n\n Simulation start time: %d/%d/%d - %d:%d:%d \n', ...
    Start_time(1),Start_time(2),Start_time(3),Start_time(4),Start_time(5),Start_time(6));
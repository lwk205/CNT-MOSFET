function CNTerrWrite(wherewri,whatis,errvalue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        CNTerrWrite.m 
% function CNTerrWrite(wherewri,whatis,errvalue)
% opens file and inital output for the carbon ananotube transistor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comb = ['\n' whatis '\n'];
fprintf(wherewri,comb);
fprintf(wherewri,'%e  ',errvalue);
fprintf(wherewri,'\n\n\n');
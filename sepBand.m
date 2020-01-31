function [En_ex,Gn_ex,Gp_ex] = sepBand(ii_band,En_old,band_old,Gn_old,Gp_old)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        sepBand.m
% function [En_ex,Gn_ex,Gp_ex] = sepBand(ii_band,En_old,band_old,Gn_old,Gp_old)
% En_ex = energies in the current band
% Gn_ex = electron polulations in the current band
% Gp_ex = hole polulations in the current band
% ii_band = current value of the band index
% En_rec = storage of energies
% band_rec = storage of band indices
% Gn_rec = storage of electron polulations
% Gp_rec = storage of hole polulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jsor = find(band_old==ii_band).';
En_ex = En_old(jsor);
Gn_ex = Gn_old(jsor,:);
Gp_ex = Gp_old(jsor,:);
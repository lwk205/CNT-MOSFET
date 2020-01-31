function [En_rec,Gn_rec,Gp_rec] = sortEn(ii_band,En_rec,band_rec,Gn_rec,Gp_rec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        sortEn.m
% function [En_rec,Gn_rec,Gp_rec] = sortEn(ii_band,En_rec,band_rec,Gn_rec,Gp_rec)
% En_rec = storage of energies
% Gn_rec = storage of electron polulations
% Gp_rec = storage of hole polulations
% ii_band = current value of the band index
% band_rec = storage of band indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jsor = find(band_rec==ii_band);
En = En_rec(jsor);
Gn = Gn_rec(jsor,:);
Gp = Gp_rec(jsor,:);
[Eso,iso] = sort(En);
Gnso = Gn(iso,:);
Gpso = Gp(iso,:);
En_rec(jsor) = Eso;
Gn_rec(jsor,:) = Gnso;
Gp_rec(jsor,:) = Gpso;

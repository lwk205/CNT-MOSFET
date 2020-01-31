% function CNTPlotResults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Title:  CNTPlotResults.m
%   function CNTPlotResults
%   Plot results produced by MOSCNT_main.m
%   The results should be already saved in 'mos_results.mat' file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
univconst                   % Define universal constants
load mos_results.mat
Ie
L = 'b-';
Ng_step = length(Vg_bias) - 1;
Nd_step = length(Vd_bias) - 1;
close all;
m = 0;

% Plot electron density
figure(1);
for ii_vg=[1:Ng_step+1-m]%1:Ng_step+1
    for ii_vd=1:Nd_step+1
        plot(XI',Ne_sub(:,ii_vg,ii_vd)*1e-2,'b-','linewidth',[2]);
        hold on;
    end
end


% grid on
h_xlabel = get(gca, 'xlabel');
h_ylabel = get(gca, 'ylabel');
title('The first subband electron density','fontsize',[18])
set(gca,'linewidth',[2],'fontsize',[12],'position',[0.15,0.2,0.74,0.7])
set(h_xlabel,'string','Positon [nm]','fontsize',[18]);
set(h_ylabel,'string','Density [cm^{-1}]','fontsize',[18]);
%hold off;


% Plot conduction/valence bands
figure(2);
for ii_vg=[1:Ng_step+1-m]%1:Ng_step+1
    for ii_vd=1:Nd_step+1
        plot(XI',E_sub(:,ii_vg,ii_vd),'b','linewidth',[2]);   hold on;
        plot(XI',E_sub(:,ii_vg,ii_vd)-2*Egh1,'r','linewidth',[2]);   hold on;
    end
end
% grid on
h_xlabel = get(gca, 'xlabel');
h_ylabel = get(gca, 'ylabel');
title('Subband energy profiles','fontsize',[18]);
set(gca,'linewidth',[2],'fontsize',[12],'position',[0.15,0.2,0.74,0.7])
set(h_xlabel,'string','Positon [nm]','fontsize',[18]);
set(h_ylabel,'string','Band edge [V]','fontsize',[18]);
%hold off;



% Plot Id-Vg
figure(3);
for ii_vd=1:Nd_step+1
    semilogy(Vg_bias(1:end)',1e6*Ie(1:end,ii_vd),L,'LineWidth',2);
    hold on
end
%axis([-inf inf 1e-8 1e2])
% grid on
h_xlabel = get(gca, 'xlabel');
h_ylabel = get(gca, 'ylabel');
set(gca,'linewidth',[2],'fontsize',[12],'position',[0.15,0.2,0.74,0.7])
set(h_ylabel,'string','I_{DS} [\muA]','fontsize',[18]);
set(h_xlabel,'string','V_{G} [V]','fontsize',[18]);
%hold off;



% Plot Id-Vd
figure(4);
for ii_vg=1:Ng_step+1
    plot(Vd_bias,1e6*Ie(ii_vg,:)','o-','LineWidth',2);
    hold on
end
% grid on
h_xlabel = get(gca, 'xlabel');
h_ylabel = get(gca, 'ylabel');
set(gca,'linewidth',[2],'fontsize',[12],'position',[0.15,0.2,0.74,0.7])
set(h_ylabel,'string','I_{DS} [\muA]','fontsize',[18]);
set(h_xlabel,'string','V_{D} [V]','fontsize',[18]);
%hold off;



% Plot Swing-Vg
figure(5);
for ii_vd=1:Nd_step+1
    S = ( Vg_bias(2:end) - Vg_bias(1:end-1) )'.*1000 ./ ( log10(Ie(2:end,ii_vd)) - log10(Ie(1:end-1,ii_vd)) );
    plot( 0.5 .* (Vg_bias(1:end-1)+Vg_bias(2:end)),abs(S),L,'LineWidth',2);
    hold on;
end
plot(0.5 .* (Vg_bias(1:end-1)+Vg_bias(2:end)),60*ones(1,length(Vg_bias)-1),'r--','LineWidth',2);
% grid on;
xlabel('V_g','fontsize',[18]);
ylabel('|Swing|, mV/dec','fontsize',[18]);
%hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plotting Gn(x,E) and LDOS(x,E)
Np = size(Gn_old,2);
for ii=1:Np
    Xcnt(ii)=acc*(3*ii/4+(1-(-1)^ii)/8-0.5);    % set atomistic grid inside coarse grid for Poisson
end
XIshift = (XI-XI(1)).';

rare = 1;
Gn_old = Gn_old(1:rare:end,:);
En_old = En_old(1:rare:end,:);


%% Plot electron distribution function
figure(5);
pcolor(Xcnt/1e-9,En_old,Gn_old);
caxis([0,0.5]);  %% Adjust upper limit to amplify the spectrum
shading interp;
colorbar;
xlabel('Position, nm','FontSize',16);
ylabel('Energy, eV','FontSize',16);
title('Gn(x,E)','FontSize',16);
hold on;
plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1),'w','linewidth',[1.5]);   hold on;
plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1)-2*Egh1,'w','linewidth',[1.5]);   hold off;


%% Plot hole distribution
Gp_old = Gp_old(1:rare:end,:);
figure(6);
pcolor(Xcnt/1e-9,En_old,Gp_old);
caxis([0,0.1]);  %% Adjust upper limit to amplify the spectrum
shading interp;
colorbar;
xlabel('Position, nm','FontSize',16);
ylabel('Energy, eV','FontSize',16);
title('Gp(x,E)','FontSize',16);
hold on;
plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1),'w','linewidth',[1]);   hold on;
plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1)-2*Egh1,'w','linewidth',[1]);   hold off;



% Plot LDOS
figure(7);
pcolor(Xcnt/1e-9,En_old,Gn_old+Gp_old);
caxis([0,0.25]);  %% Adjust upper limit to amplify the spectrum
shading interp;
colorbar;
xlabel('Position, nm','FontSize',16);
ylabel('Energy, eV','FontSize',16);
title('LDOS(x,E)','FontSize',16);
hold on;
plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1),'w','linewidth',[1]);   hold on;
plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1)-2*Egh1,'w','linewidth',[1]);   hold off;



%% works only for: (WKB_flag = 'rigor' && scat_flag = 'negf') combination
%%% plot I(E) spectrum
% rare = 1;
% myEE = myEE(1:rare:end);
% myTr = myTr(1:rare:end);
% Trraspr = Trraspr(1:rare:end,:);
% [myEE_new my_index] = sort(myEE);
% myTr_new = myTr(my_index);
% Trraspr_new = Trraspr(my_index,:);
% 
% % plot I(E,x) spectrum
% minI = min(min(Trraspr_new))
% maxI = max(max(Trraspr_new))
% novT = real(log10(Trraspr_new));
% jneg = find(Trraspr_new<0);
% %novT(jneg) = -novT(jneg);
% 
% figure(9);
% pcolor(Xcnt(1:end-1)/1e-9,myEE_new(1:rare:end).',novT(1:rare:end,:));
% caxis([-8,0]);  %% Adjust limits to amplify the spectrum
% shading interp;
% colorbar;
% xlabel('Position, nm','FontSize',16);
% ylabel('Energy, eV','FontSize',16);
% title('Log_{10} I(x,E)','FontSize',16);
% hold on;
% plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1),'w','linewidth',[1.5]);   hold on;
% plot(XIshift,E_sub(:,Ng_step+1,Nd_step+1)-2*Egh1,'w','linewidth',[1.5]);
% hold off;

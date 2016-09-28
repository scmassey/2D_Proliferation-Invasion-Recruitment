
% run_pir_sim(eta_c,Dwg_factor,dt,dayspersave,t_finish,...
%     radiussave,radius_final,savefile,center_ind1,...
%     center_ind2)

clear all 
close all
clc

% eta_c - defined in loop below
% Dwg_factor - defined in loop below
species='mouse';
start = 'gray';
% dt - defined in loop below (to handle stiffness)
dayspersave=10; %set 0 if save by radius; was 10;% not actually days, but time steps - since dt is 0.5, this is 5 days
% t_finish- defined below to match with dt; =100; % 100 steps = 50 days; 20 days = 40 steps
radiussave=0.05;% every half mm
radius_final=0.5;% stop at half cm (since brain radius is 1 cm, above this seems too high)
% savefile - defined in loop below

center_ind1 = 87; % starts in gray;  % center_ind1 = 98; % starts in white
center_ind2 = 117;                   % center_ind2 = 110;

ivals = [0, 10^(-6),10^(-5),10^(-4)]; 
jvals = [5,10,50,100];
for i=1:length(ivals)
    for j=1:length(jvals)
        eta_c=ivals(i);
        Dwg_factor=jvals(j);
        if i==4 && j==3
            dt=0.25; % days
            t_finish = 120;
        elseif i==4 && j==4
            dt=0.2; % days
            t_finish = 150; 
        else 
            dt=0.5; % days
            t_finish = 60;
        end;
        savestr=['2dPIR_sim_SeptFIXED_',species,'_',start,'_',num2str(center_ind1),'_',num2str(center_ind2),'_eta_',num2str(i),'_DwDg_',num2str(jvals(j)),'_times_newR0'];
        run_pir_sim_newR0(eta_c,Dwg_factor,species,dt,dayspersave,t_finish,radiussave,...
            radius_final,savestr,center_ind1,center_ind2)
    end;
end;
% To visualize results, use plot_2dpir_sims_set.m
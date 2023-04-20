%% get_hybrid_model_results.m:
% Script to generate results from the VOC importation model
% and hybrid stochastic-parsimonious transmission model
%--------------------------------------------------------------------------
clear

%% Generate mex files for run_gillespie_model_for_mex
R_eff_wildtype = 3.51;
threshold_prevalence = 100;
relsuscrec = 0;  % If infected by UK resident variants previously, have immunity to VOC
relsuscvacc = [1 0.35 0.25 0]; %Susceptibilties for [No vacc, AZ, Pfizer, VOC targeted];
time_horizon = 365; % Set duration to run the outbreak for
effective_imports = 1;
VOC_rel_trans = 1;
codegen run_gillespie_model_for_mex -args {R_eff_wildtype,VOC_rel_trans,relsuscvacc,relsuscrec,threshold_prevalence,time_horizon,effective_imports}

%--------------------------------------------------------------------------
%% Set up runs varying four things: imports, immune escape, transmissibility and R_eff_wildtype
%--------------------------------------------------------------------------
% Set global parameters
threshold_prevalence = 100;
relsuscrec = 0;  % If infected by wildtype previously, have immunity to VOC
relsuscvacc = [1 0.35 0.25 0]; %Susceptibilties for [No vacc, AZ, Pfizer, VOC targeted];
time_horizon = 365; % Set duration to run the outbreak for

effective_imports_over = 0.02:0.02:0.4;
VOC_rel_trans_over =  0.5:0.1:1.5;
R_excl_immun_wildtype_over = [3];   % R excluding immunity for the resident variants
relative_suscept_over = 0.75;

% Set storage arrays
epidemic_prob = zeros(numel(effective_imports_over),numel(VOC_rel_trans_over),numel(R_excl_immun_wildtype_over),numel(relative_suscept_over));
reach_thresh_time = zeros(numel(effective_imports_over),numel(VOC_rel_trans_over),numel(R_excl_immun_wildtype_over),numel(relative_suscept_over));
all_output_at_thresh = zeros(numel(effective_imports_over),numel(VOC_rel_trans_over),numel(R_excl_immun_wildtype_over),numel(relative_suscept_over),1000,17);

% Run the scenarios
length_effective_imports = numel(effective_imports_over);
length_VOC_rel_trans = numel(VOC_rel_trans_over);
length_R_excl_immun_wildtype = numel(R_excl_immun_wildtype_over);
length_relative_suscept_over = numel(relative_suscept_over);
parfor effective_imports_itr = 1:length_effective_imports
    disp(['it = ',mat2str(effective_imports_itr),' out of ',mat2str(length_effective_imports)])
    effective_imports = effective_imports_over(effective_imports_itr);
    
    for VOC_rel_trans_itr = 1:length_VOC_rel_trans
        VOC_rel_trans = VOC_rel_trans_over(VOC_rel_trans_itr);
        display(['itr = ',mat2str(VOC_rel_trans_itr),' out of ',mat2str(length_VOC_rel_trans)])
        
        for R_excl_immun_wildtype_itr = 1:length_R_excl_immun_wildtype
            R_excl_immun_wildtype = R_excl_immun_wildtype_over(R_excl_immun_wildtype_itr);
            
            
            for relative_suscept_over_itr = 1:length_relative_suscept_over
                relative_suscept = relative_suscept_over(relative_suscept_over_itr);
                
                % Pass parameters to Gillespie stochastic simulation
                [epidemic_prob(effective_imports_itr,VOC_rel_trans_itr,R_excl_immun_wildtype_itr,relative_suscept_over_itr),...
                    reach_thresh_time(effective_imports_itr,VOC_rel_trans_itr,R_excl_immun_wildtype_itr,relative_suscept_over_itr),...
                    all_output_at_thresh(effective_imports_itr,VOC_rel_trans_itr,R_excl_immun_wildtype_itr,relative_suscept_over_itr,:,:)] = run_gillespie_model_for_mex_mex(R_excl_immun_wildtype,...
                    VOC_rel_trans,...
                    (1-(1-relsuscvacc)*relative_suscept),...
                    (1-relative_suscept),...
                    threshold_prevalence,...
                    time_horizon,...
                    effective_imports);
            end
        end
    end
end

% % Save the outputs
% save('gillespie_model_outputs/gillespie_model_varying_all4.mat','epidemic_prob','reach_thresh_time','all_output_at_thresh',...
%     'effective_imports_over','VOC_rel_trans_over','R_excl_immun_wildtype_over','relative_suscept_over')

%%
%--------------------------------------------------------------------------
%% Parsimonious model runs using stochastic initial conditions
%--------------------------------------------------------------------------
% Load initial conditions from file
load('MAT_files/gillespie_model_varying_all4.mat')

% Decide what to change
Tidx = 1; Eidx = 2:9; Iidx = 10:17; % Indices of time columns, E columns (8 TaB), and I columns (8 TaB)

% Set up loop variables
new_vaccine_intro_date_over = [datenum(2021,6,1) datenum(2021,7,1) datenum(2021,8,1) datenum(2021,9,1) datenum(2021,10,1) datenum(2021,11,1)];
ordering = {'effective_imports_over','VOC_rel_trans_over','R_excl_immun_wildtype_over','relative_suscept_over'};
    % 'new_vaccine_intro_date_over' is also possible
plot_over_x_name = 'effective_imports_over';
plot_over_y_name = 'relative_suscept_over';
plot_over_z_name = 'new_vaccine_intro_date_over';
plot_over_x = eval(plot_over_x_name);
plot_over_y = eval(plot_over_y_name);
plot_over_z = eval(plot_over_z_name);

% set default values for other parameters
effective_imports_pos = find(effective_imports_over==0.2); % set effective imports per day
VOC_rel_trans_pos = find(VOC_rel_trans_over==1);
R_excl_immun_wildtype_pos = find(R_excl_immun_wildtype_over==3);
relative_suscept_pos = find(relative_suscept_over==1);

% find start date
parameters = make_parameters();
start_date = parameters.date1;
maxT = parameters.maxT;

epidemic_size_gillespie = NaN(length(plot_over_x),length(plot_over_y),length(plot_over_z),1000,3);
peak_height_gillespie = NaN(length(plot_over_x),length(plot_over_y),length(plot_over_z),1000,3);
peak_time_gillespie = NaT(length(plot_over_x),length(plot_over_y),length(plot_over_z),1000,3);
for iterate_x = 1:length(plot_over_x)
    clear changed_parameters
    display(['it = ',mat2str(iterate_x),' out of ',mat2str(length(plot_over_x))])
    
    % Set up changed_parameters structure
    changed_parameters = set_parameters(plot_over_x_name,plot_over_x,iterate_x,struct());
    
    % Specify that we will be specifying initial condition distribution
    changed_parameters.specify_distribution = true;
    
    % Iterate over immune escape parameters
    for iterate_y = 1:length(plot_over_y)
        changed_parameters = set_parameters(plot_over_y_name,plot_over_y,iterate_y,changed_parameters);
        for iterate_z = 1:length(plot_over_z)
            changed_parameters = set_parameters(plot_over_z_name,plot_over_z,iterate_z,changed_parameters);
            % find initial VOC
            Index = {effective_imports_pos,VOC_rel_trans_pos,R_excl_immun_wildtype_pos,relative_suscept_pos,':',':'};
            
            % find matrix for plot_over_x and y
            if sum(strcmp(ordering,plot_over_x_name))>0
                Index{strcmp(ordering,plot_over_x_name)} = iterate_x;
            end
            if sum(strcmp(ordering,plot_over_y_name))>0
                Index{strcmp(ordering,plot_over_y_name)} = iterate_y;
            end
            
            % Use indexes to look up relevant initial conditions from
            % Gillespie simulation outputs
            vals = squeeze(all_output_at_thresh(Index{:}))';
            
            % For those simulations that achieved the specified threshold
            % value (Plim), get the times the threshold value was attained.
            % Column 1 of vals corresponds to time to hit Plim.
            times = vals(1,~isnan(vals(1,:)));
            
            for intro_itr = 1:length(times)
                if times(intro_itr)<(maxT-1)
                    changed_parameters.VOC_imp_date = start_date + times(intro_itr);
                    changed_parameters.VOC_imp_distribution(1,1,:) = vals(Eidx(1:4),intro_itr)/56e6;
                    changed_parameters.VOC_imp_distribution(2,1,:) = vals(Eidx(5:8),intro_itr)/56e6;
                    changed_parameters.VOC_imp_distribution(1,2,:) = vals(Iidx(1:4),intro_itr)/56e6;
                    changed_parameters.VOC_imp_distribution(2,2,:) = vals(Iidx(5:8),intro_itr)/56e6;
                    parameters = make_parameters(changed_parameters);
                    
                    % run the simple model
                    [t,pop_out_gillespie,parameters_out_gillespie,outputs_gillespie] = run_simple_vaccines_mex(parameters);
                    
                    % find outputs
                    [epidemic_size_gillespie(iterate_x,iterate_y,iterate_z,intro_itr,:),...
                        peak_height_gillespie(iterate_x,iterate_y,iterate_z,intro_itr,:),...
                        peak_time_gillespie(iterate_x,iterate_y,iterate_z,intro_itr,:)] = process_outputs(parameters_out_gillespie,outputs_gillespie,pop_out_gillespie);
                end
            end
        end
    end
end
% save('VOC_E_runs.mat','epidemic_size_gillespie',...
%                         'peak_height_gillespie',...
%                          'peak_time_gillespie',...
%                          'parameters',...
%                          'plot_over_x',...
%                          'plot_over_x_name',...
%                          'plot_over_y',...
%                          'plot_over_y_name',...
%                          'plot_over_z',...
%                          'plot_over_z_name')

%--------------------------------------------------------------------------
%% SUPPORTING FUNCTIONS %%
%--------------------------------------------------------------------------                                   

% Compute epidemic size, peak size and timing of the peak.
function [epidemic_size,peak_height,peak_time] = process_outputs(parameters,outputs,pop_out)
    
    % Get field names in parameter structure
    % Assign each field entry to its own separate variable
    names = fieldnames(parameters);
    for i=1:length(names)
        eval([cell2mat(names(i)),' = parameters.',cell2mat(names(i)),';']);
    end
    
    % Get field names in outputs structure
    % Assign each field entry to its own separate variable
    names = fieldnames(outputs);
    for i=1:length(names)
        eval([cell2mat(names(i)),' = outputs.',cell2mat(names(i)),';']);
    end
    
    % Get peaks in prevalence pk_H_XXX (and the timing, pk_t_XXX)
    % for resident variant, VOC & both variants summed
    [pk_H_UK,pk_t_UK] = max(I_UK);
    [pk_H_VOC,pk_t_VOC] = max(I_VOC);
    [pk_H_both,pk_t_both] = max(I_UK+I_VOC);
    
    % Assign prevalence summary statistics to output variables
    peak_time = [dates(pk_t_UK),dates(pk_t_VOC),dates(pk_t_both)];
    peak_height = [pk_H_UK,pk_H_VOC,pk_H_both];
    epidemic_size = [sum(pop_out(4,:,:,end),'all'),sum(pop_out(:,4,:,end),'all'),sum(pop_out(4,:,:,end),'all')+sum(pop_out(:,4,:,end),'all')];
end

%--------------------------------------------------------------------------
%% Supporting function for iterating over arbitrary variables
%--------------------------------------------------------------------------
function changed_parameters = set_parameters(plot_over_x_name,plot_over_x,iterate_x,changed_parameters)
if strcmp(plot_over_x_name,'VOC_rel_trans_over')
    changed_parameters.VOC_vs_UK = plot_over_x(iterate_x);
elseif strcmp(plot_over_x_name,'relative_suscept_over')
    % - Vaccine efficacy scaling
    changed_parameters.e_aVOC_scaling =  plot_over_x(iterate_x);
    changed_parameters.e_pVOC_scaling =  plot_over_x(iterate_x);
    
    % - Cross-immunity
    changed_parameters.s_VOC = 1- plot_over_x(iterate_x); % susceptibility to VOC for resident variant recovereds
    changed_parameters.s_UK = 1- plot_over_x(iterate_x); % susceptibility to resident variants for VOC recovereds
elseif strcmp(plot_over_x_name,'new_vaccine_intro_date_over')
    parameters = make_parameters();
    changed_parameters.vaccine_changeover_week = floor((plot_over_x(iterate_x) - parameters.date1)/7);
    changed_parameters.prioritise_unvaccinated = 4; % prioritise unvaccinated, then either AZ or P 
    % changed_parameters.prioritise_unvaccinated = 5; % prioritise vaccinated, then unvaccinated group. 
end
end
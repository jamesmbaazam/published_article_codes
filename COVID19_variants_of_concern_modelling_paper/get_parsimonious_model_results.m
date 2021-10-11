%% get_parsimonious_model_results.m:
% Script to generate outputs that were used to produce the
% figures associated with the parsimonious SARS-CoV-2
% transmission model and the stochastic importation model.
%--------------------------------------------------------------------------

clear

%% Set global flag variables
make_mex_flag = true;

%% Make mex: run this the first time to make the mex file
if make_mex_flag == true
    %%
    clear changed_parameters
    changed_parameters.VOC_imp_size = 0;
    parameters = make_parameters(changed_parameters);
    codegen run_simple_vaccines -args {parameters}
end

%--------------------------------------------------------------------------
%% NO VOC RUN %%
%--------------------------------------------------------------------------
clear changed_parameters
changed_parameters.VOC_imp_size = 0;
parameters = make_parameters(changed_parameters);
[t,no_VOC_pop_out,no_VOC_parameters,no_VOC_outputs] = run_simple_vaccines_mex(parameters);

% Get effective R for each variant
[~,ReffUK_no_VOC] = get_Reff(no_VOC_pop_out,no_VOC_parameters);

% Get infectious temporal profiles from outputs data
I_UK_no_VOC = no_VOC_outputs.I_UK;

% Get vaccination coverage data
vacc_coverage_data = squeeze(sum(no_VOC_pop_out(:,:,2:4,:),[1,2,3]))*100;

%--------------------------------------------------------------------------
%% BASE PARAMETER RUN %%
%--------------------------------------------------------------------------

%% Get base level data
parameters = make_parameters();
[t,base_pop_out,base_parameters,base_outputs] = run_simple_vaccines_mex(parameters);

%--------------------------------------------------------------------------
%% VOC SCENARIO RUNS.
%--------------------------------------------------------------------------

% Set up VOC parameters
VOC_vs_UK_varies = [1.5 1   0.8 1   1]; % Relative transmissibility of VOC
efficacy_varies =  [1   0.75 0.75 0.75 1 ]; % Vaccine efficacy
s_varies =         [1   0.75 0.75 1   0.75]; % Natural immunity efficacy
n_VOCs = numel(VOC_vs_UK_varies);

% Set up introduction dates for VOCs
n_intro_dates = datenum(2021,11,1) - datenum(2021,5,17) + 1;

% Set up storage arrays for effective R
ReffVOC_default_runs = zeros(n_VOCs,parameters.maxT+1,n_intro_dates);
ReffUK_default_runs= zeros(n_VOCs,parameters.maxT+1,n_intro_dates);

% Set up storage arrays for infectious temporal profiles
I_VOC_default_runs = zeros(parameters.maxT+1,n_VOCs,n_intro_dates);
I_UK_default_runs = zeros(parameters.maxT+1,n_VOCs,n_intro_dates);

% Initialise output variables for outbreak size and peak in infection
epidemic_size_default_runs = zeros(n_intro_dates,n_VOCs,3);
peak_height_default_runs = zeros(n_intro_dates,n_VOCs,3);

% Iterate over introduction dates for VOC
for jj = 1:n_intro_dates
    clear changed_parameters
    % changed_parameters.VOC_imp_date = VOC_intro_dates(jj);
    changed_parameters.VOC_imp_date = datenum(2021,5,17) + (jj-1);
    
    % Iterate over each VOC
    % For each, compute effective R and epidemiological outputs
    for ii = 1:n_VOCs

        % Set up VOC parameters
        % - Relative transmissibility
        changed_parameters.VOC_vs_UK = VOC_vs_UK_varies(ii);

        % - Vaccine efficacy scaling
        changed_parameters.e_aVOC_scaling = efficacy_varies(ii);
        changed_parameters.e_pVOC_scaling = efficacy_varies(ii);

        % - Cross-immunity
        changed_parameters.s_VOC = 1-s_varies(ii); % susceptibility to VOC for resident variants recovereds
        changed_parameters.s_UK = 1-s_varies(ii); % susceptibility to resident variants for VOC recovereds

        % Run the model
        parameters = make_parameters(changed_parameters);
        [t,pop_out_default_runs,parameters_default_runs,outputs_default_runs] = run_simple_vaccines_mex(parameters);
        [epidemic_size_default_runs(jj,ii,:),peak_height_default_runs(jj,ii,:),peak_time_default_runs(jj,ii,:)] = process_outputs(parameters_default_runs,outputs_default_runs,pop_out_default_runs);        

        % Get effective R for each variant
        [ReffVOC_default_runs(ii,:,jj),ReffUK_default_runs(ii,:,jj)] = get_Reff(pop_out_default_runs,parameters_default_runs);

        % Get infectious temporal profiles from outputs data
        I_VOC_default_runs(:,ii,jj) = outputs_default_runs.I_VOC;
        I_UK_default_runs(:,ii,jj) = outputs_default_runs.I_UK;
    end
end  

%--------------------------------------------------------------------------
%% VOC SCENARIO RUNS. R with immunity
%% Figures 1(c), S2(b)
%--------------------------------------------------------------------------                       
clear changed_parameters

% Get effective R for each VOC and the resident variants
for ii=1:n_VOCs
    clear changed_parameters    
    
    % no VOC    
    changed_parameters.VOC_imp_size = 0;
    
    % Set up VOC parameters
    % - Relative transmissibility
    changed_parameters.VOC_vs_UK = VOC_vs_UK_varies(ii);
    
    % - Vaccine efficacy scaling
    changed_parameters.e_aVOC_scaling = efficacy_varies(ii);
    changed_parameters.e_pVOC_scaling = efficacy_varies(ii);
    
    % - Cross-immunity
    changed_parameters.s_VOC = 1-s_varies(ii); % susceptibility to VOC for resident variants recovereds
    changed_parameters.s_UK = 1-s_varies(ii); % susceptibility to resident variants for VOC recovereds
    
    parameters = make_parameters(changed_parameters);
    [t,no_VOC_pop_out,no_VOC_parameters,no_VOC_outputs] = run_simple_vaccines_mex(parameters);
    
    % Get R with immunity values
    [R0VOC{ii},R0UK] = get_Reff(no_VOC_pop_out,no_VOC_parameters);
end

%--------------------------------------------------------------------------
%% SENSITIVITY HEATMAPS (FIGURES 1B & 1D) %%
%--------------------------------------------------------------------------

%% Change relative VOC_vs_UK transmissibility and AZ vaccine efficacy against VOC together

% Specify values to test
% efficacy_varies = 0.9:-0.05:0.1;
VOC_efficacy_scaling = 0.5:0.05:1;
VOC_vs_UK_varies = 0.5:0.1:1.5;

% Initialise output variables
epidemic_size_1 = zeros(length(VOC_efficacy_scaling),length(VOC_vs_UK_varies),3);
peak_height_1 = zeros(length(VOC_efficacy_scaling),length(VOC_vs_UK_varies),3);

% Iterate over each parameter combination
for ii=1:length(VOC_efficacy_scaling)
    clear changed_parameters
    disp(['ii = ',mat2str(ii),', out of ',mat2str(length(VOC_efficacy_scaling))])
    
    % Update vaccine efficacy scaling
    changed_parameters.e_aVOC_scaling = VOC_efficacy_scaling(ii);
    changed_parameters.e_pVOC_scaling = VOC_efficacy_scaling(ii);

    % Update susceptibility given prior infection
    % Equivalent to (1 - efficacy)
    changed_parameters.s_UK = 1 - VOC_efficacy_scaling(ii);
    changed_parameters.s_VOC = 1 - VOC_efficacy_scaling(ii);

    % Perform sweep over relative transmissibility values
    for jj=1:length(VOC_vs_UK_varies)
        changed_parameters.VOC_vs_UK = VOC_vs_UK_varies(jj);
        parameters = make_parameters(changed_parameters);
        [t,pop_out_1{ii,jj},parameters_1{ii,jj},outputs_1{ii,jj}] = run_simple_vaccines_mex(parameters);
        [epidemic_size_1(ii,jj,:),peak_height_1(ii,jj,:),peak_time_1(ii,jj,:)] = process_outputs(parameters_1{ii,jj},outputs_1{ii,jj},pop_out_1{ii,jj});
    end
end

%--------------------------------------------------------------------------
%% R with immunity for VOC introduced on different calendar dates and with
%% various transmissibility & immune escape profiles
%% Used for Figure S5
%--------------------------------------------------------------------------

% Made as snapshots in time for different calendar dates.
% Date for seeding initial VOC infecteds & time snapshots to compute
% relative R
VOC_imp_date_varies = [datenum(2021,5,17) datenum(2021,8,1) datenum(2021,11,1)];
t_snapshots = [datenum(2021,5,17) datenum(2021,8,1) datenum(2021,11,1)];

% Transmissibility varied
VOC_vs_UK_varies = 0.5:0.1:1.5;

% Vaccine efficacy varied
VOC_efficacy_scaling = 0.5:0.05:1;

% Initialise output variables
rel_eff_R_heatmap = zeros(length(VOC_efficacy_scaling),length(VOC_vs_UK_varies),numel(t_snapshots));
VOC_eff_R_heatmap = zeros(length(VOC_efficacy_scaling),length(VOC_vs_UK_varies),numel(t_snapshots));
resident_variants_eff_R_heatmap = zeros(length(VOC_efficacy_scaling),length(VOC_vs_UK_varies),numel(t_snapshots));

% Iterate over each parameter combination
% For specified timepoints/seeding of initial VOC infecteds, get the relative R
for kk = 1:numel(t_snapshots)
    clear changed_parameters
    disp(['kk = ',mat2str(kk),', out of ',mat2str(numel(t_snapshots))])
    changed_parameters.VOC_imp_date = VOC_imp_date_varies(kk);
    t_snapshot_slice_idx = t_snapshots(kk) - parameters.date1 + 1;
    for ii=1:length(VOC_efficacy_scaling)
        % - Vaccine efficacy scaling
        changed_parameters.e_aVOC_scaling = VOC_efficacy_scaling(ii);
        changed_parameters.e_pVOC_scaling = VOC_efficacy_scaling(ii);

        % - Cross-immunity
        changed_parameters.s_VOC = 1-VOC_efficacy_scaling(ii); % susceptibility to VOC for resident variants recovereds
        changed_parameters.s_UK = 1-VOC_efficacy_scaling(ii); % susceptibility to resident variants for VOC recovereds

        for jj=1:length(VOC_vs_UK_varies)
            % Run the model for the current parameter set
            changed_parameters.VOC_vs_UK = VOC_vs_UK_varies(jj);
            parameters = make_parameters(changed_parameters);
            [t,pop_out,parameters,outputs] = run_simple_vaccines_mex(parameters);
               
            % Get effective R for each variant
            % Compute relative values
            [ReffVOC,ReffUK] = get_Reff(pop_out,parameters);
            rel_eff_R_heatmap(ii,jj,kk) =  ReffVOC(t_snapshot_slice_idx)/ReffUK(t_snapshot_slice_idx);
            VOC_eff_R_heatmap(ii,jj,kk) =  ReffVOC(t_snapshot_slice_idx);
            resident_variants_eff_R_heatmap(ii,jj,kk) =  ReffUK(t_snapshot_slice_idx);
        end
    end
end

% Get threshold measure of relative effective R > 1
threshold_rel_eff_R_heatmap = 0*rel_eff_R_heatmap;
threshold_rel_eff_R_heatmap(rel_eff_R_heatmap>1 & VOC_eff_R_heatmap>1) = 1;

% Store outputs in cell for use in plotting
effective_R_data = {rel_eff_R_heatmap,VOC_eff_R_heatmap,resident_variants_eff_R_heatmap,threshold_rel_eff_R_heatmap};

%--------------------------------------------------------------------------
%% TRANSMISSION BLOCKING SENSITIVITY RUNS  %%
%% Used in Figures S6 & S7
%--------------------------------------------------------------------------

% Options for varying the transmission blocking assumptions
transmission_blocking_val = [0 0.25 0.5];
n_transmission_blocking_vals = numel(transmission_blocking_val);

% Set up VOC parameters
VOC_vs_UK_varies = [1.5 1   0.8]; % Relative transmissibility of VOC
efficacy_varies =  [1   0.75 0.75]; % Vaccine efficacy
s_varies =         [1   0.75 0.75]; % Natural immunity efficacy
n_VOCs = numel(VOC_vs_UK_varies);

% Set up storage arrays 
ReffVOC_trans_block = zeros(n_transmission_blocking_vals,n_VOCs,parameters.maxT+1);
ReffUK_trans_block = zeros(n_transmission_blocking_vals,n_VOCs,parameters.maxT+1);
ReffUK_no_VOC_trans_block = zeros(n_transmission_blocking_vals,parameters.maxT+1);

% Set up storage arrays for infectious temporal profiles
I_VOC_trans_block = zeros(parameters.maxT+1,n_VOCs,n_transmission_blocking_vals);
I_UK_trans_block = zeros(parameters.maxT+1,n_VOCs,n_transmission_blocking_vals);
I_UK_no_VOC_trans_block = zeros(parameters.maxT+1,n_transmission_blocking_vals);

% Iterate over transmission blocking values
% For each, get effective R for the VOC
for jj = 1:n_transmission_blocking_vals
    disp(['jj = ',mat2str(jj),', out of ',mat2str(numel(transmission_blocking_val))])
    clear changed_parameters
    propn_remaining_transmission = 1 - transmission_blocking_val(jj);
    changed_parameters.propn_transmission_a = propn_remaining_transmission;
    changed_parameters.propn_transmission_p = propn_remaining_transmission;
    changed_parameters.propn_transmission_n = propn_remaining_transmission;
    changed_parameters.propn_transmission_priorinf = propn_remaining_transmission;
    for ii = 1:n_VOCs
        % Set up VOC parameters
        changed_parameters.VOC_vs_UK = VOC_vs_UK_varies(ii);
        changed_parameters.e_aVOC_scaling = efficacy_varies(ii);
        changed_parameters.e_pVOC_scaling = efficacy_varies(ii);
        changed_parameters.s_VOC = 1-s_varies(ii); % susceptibility to VOC variant for UK recovereds
        changed_parameters.s_UK = 1-s_varies(ii); % susceptibility to VOC variant for UK recovereds

        % Run the model
        parameters = make_parameters(changed_parameters);
        [t,pop_out_trans_block,parameters_trans_block,outputs_trans_block] = run_simple_vaccines_mex(parameters);
        
        % Get effective R for each variant
        % Compute relative values
        [ReffVOC_trans_block(jj,ii,:),ReffUK_trans_block(jj,ii,:)] = get_Reff(pop_out_trans_block,parameters_trans_block);
        
        % Get infectious temporal profiles from outputs data
        I_VOC_trans_block(:,ii,jj) = outputs_trans_block.I_VOC;
        I_UK_trans_block(:,ii,jj) = outputs_trans_block.I_UK;
    end
    
    % Run with no VOCs
    changed_parameters.VOC_imp_size = 0;
    parameters = make_parameters(changed_parameters);
    [t,no_VOC_trans_block_pop_out,no_VOC_trans_block_parameters,no_VOC_trans_block_outputs] = run_simple_vaccines_mex(parameters);

    % Get effective R for each variant
    [~,ReffUK_no_VOC_trans_block(jj,:)] = get_Reff(no_VOC_trans_block_pop_out,no_VOC_trans_block_parameters);

    % Get infectious temporal profiles from outputs data
    I_UK_no_VOC_trans_block(:,jj) = no_VOC_trans_block_outputs.I_UK;
end

%--------------------------------------------------------------------------
%% SAVE OUTPUTS TO FILE %%
%--------------------------------------------------------------------------
save('MAT_files/parsimonious_model_results.mat',...
                                       'vacc_coverage_data',...
                                       'base_parameters',... 
                                       'base_outputs',...
                                       'I_UK_no_VOC',...
                                       'no_VOC_outputs',...
                                       'I_VOC_default_runs',...
                                       'I_UK_default_runs',...
                                       'outputs_default_runs',...
                                       'ReffVOC_default_runs',...
                                       'ReffUK_default_runs',...
                                       'R0VOC',...
                                       'R0UK',...
                                       'epidemic_size_1',...
                                       'peak_height_1',...
                                       'peak_time_1',...
                                       'epidemic_size_default_runs',...
                                       'peak_height_default_runs',...
                                       'peak_time_default_runs',...
                                       'I_VOC_trans_block',...
                                       'I_UK_trans_block',...
                                       'I_UK_no_VOC_trans_block',...
                                       'ReffVOC_trans_block',...
                                       'ReffUK_trans_block',...
                                       'ReffUK_no_VOC_trans_block',...
                                       'outputs_trans_block',...
                                       'effective_R_data');
                                                         
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

% Compute effective R (R with immunity)
function [ReffVOC,ReffUK] = get_Reff(pop_out,parameters)

    % Get number of timesteps
    t = 1:size(pop_out,4);

    %% calculate immunity to the VOC variant
    clear VOC_perc_immun_pop_out
    VOC_perc_immun_pop_out = 0*pop_out;

    % Potentially have some VOC immunity if
    % - Currently infected with VOC or a recovered from VOC infection
    % - Susceptible to VOC & currently infected with resident variants (no co-infection allowed)
    % - Vaccinated (no prior infection with any variant)
    % - Unvaccinated, not previously infected by VOC, previously infected resident variants
    % - Vaccinated, not previously infected by VOC, previously infected by resident variants

    % if you're VOC recovered or currently infected:
    VOC_perc_immun_pop_out(:,2:4,:,:) = pop_out(:,2:4,:,:);

    % Susceptible to resident variants & currently infected with VOC 
    % (no co-infection allowed)
    VOC_perc_immun_pop_out(2:3,1,:,:) = pop_out(2:3,1,:,:);

    % Vaccianted with no prior infection
    VOC_perc_immun_pop_out(1,1,2,:) = pop_out(1,1,2,:)*(1-parameters.e_aVOC); % AZ;
    VOC_perc_immun_pop_out(1,1,3,:) = pop_out(1,1,3,:)*(1-parameters.e_pVOC); % Pfizer
    VOC_perc_immun_pop_out(1,1,4,:) = pop_out(1,1,4,:)*(1-parameters.e_nVOC); % New vaccine

    % if you're UK recovered and unvaccinated & susceptible to VOC variant:
    VOC_perc_immun_pop_out(4,1,1,:) = pop_out(4,1,1,:)*(1-parameters.s_VOC);

    % if you're both UK recovered and vaccinated:
    VOC_perc_immun_pop_out(4,1,2,:) = pop_out(4,1,2,:)*(1-min([parameters.s_VOC,parameters.e_aVOC]));
    VOC_perc_immun_pop_out(4,1,3,:) = pop_out(4,1,3,:)*(1-min([parameters.s_VOC,parameters.e_pVOC]));
    VOC_perc_immun_pop_out(4,1,4,:) = pop_out(4,1,4,:)*(1-min([parameters.s_VOC,parameters.e_nVOC]));

    % Get percentage immune at each timestep (sum over dimension 1 through 3, with dimension 4 having a slice per timestep)
    VOC_perc_immun = sum(VOC_perc_immun_pop_out,[1,2,3]);

    %% calculate immunity to the UK variant
    clear UK_perc_immun_pop_out
    UK_perc_immun_pop_out = 0*pop_out;

    % Potentially have some resident variants immunity if
    % - Currently infected with resident variants or a recovered from resident variants infection
    % - Susceptible to resident variants & currently infected with VOC (no co-infection allowed)
    % - Vaccinated (no prior infection with any variant)
    % - Unvaccinated, not previously infected by resident variants, previously infected by VOC
    % - Vaccinated, not previously infected by resident variants, previously infected by VOC

    % if you're UK recovered or currently infected:
    UK_perc_immun_pop_out(2:4,:,:,:) = pop_out(2:4,:,:,:);

    % Susceptible to resident variants & currently infected with VOC 
    % (no co-infection allowed)
    UK_perc_immun_pop_out(1,2:3,:,:) = pop_out(1,2:3,:,:);

    % Vaccianted with no prior infection
    UK_perc_immun_pop_out(1,1,2,:) = pop_out(1,1,2,:)*(1-parameters.e_aUK); % AZ;
    UK_perc_immun_pop_out(1,1,3,:) = pop_out(1,1,3,:)*(1-parameters.e_pUK); % Pfizer
    UK_perc_immun_pop_out(1,1,4,:) = pop_out(1,1,4,:)*(1-parameters.e_nUK); % New vaccine

    % if you're VOC recovered and unvaccinated & susceptible to UK variant:
    UK_perc_immun_pop_out(1,4,1,:) = pop_out(1,4,1,:)*(1-parameters.s_UK);

    % if you're both VOC recovered and vaccinated:
    UK_perc_immun_pop_out(1,4,2,:) = pop_out(1,4,2,:)*(1-min([parameters.s_UK,parameters.e_aUK]));
    UK_perc_immun_pop_out(1,4,3,:) = pop_out(1,4,3,:)*(1-min([parameters.s_UK,parameters.e_pUK]));
    UK_perc_immun_pop_out(1,4,4,:) = pop_out(1,4,4,:)*(1-min([parameters.s_UK,parameters.e_nUK]));

    % Get percentage immune at each timestep (sum over dimension 1 through 3, with dimension 4 having a slice per timestep)
    UK_perc_immun = sum(UK_perc_immun_pop_out,[1,2,3]);
    temp = squeeze(UK_perc_immun);
    %temp(end)
    %% Get an "average" infectious individual

    % Get total number of infecteds at each time point
    I_VOC_per_t = squeeze(sum(sum(pop_out(:,3,:,:),1),3));
    I_UK_per_t = squeeze(sum(sum(pop_out(3,:,:,:),2),3));

    % Get force of infection, weighted by any transmission blocking effect
    [weighted_I_UK,weighted_I_VOC] = compute_weighted_force_of_infection(pop_out,parameters);

    % Get scaling factor for transmissibility. Domain [0,1] 
    %   1 corresponds to unmodified, no transmission blocking.
    %   0 corresponds to all transmission blocked.
    avg_trans_VOC_inf = weighted_I_VOC./I_VOC_per_t;
    avg_trans_VOC_inf(I_VOC_per_t==0) = 1; %Prior to VOC being seeded, reset to 1
    avg_trans_UK_inf = weighted_I_UK./I_UK_per_t;

    %% Compute R effective
    bet_VOC = ones(size(pop_out,4),1)*parameters.beta_VOC_changes(end);
    bet_UK = ones(size(pop_out,4),1)*parameters.beta_UK_changes(end);
    for i=length(parameters.change_days):-1:1
        bet_VOC(t<=parameters.change_days(i))=parameters.beta_VOC_changes(i);
        bet_UK(t<=parameters.change_days(i))=parameters.beta_UK_changes(i);
    end
    ReffVOC = (bet_VOC./parameters.gam).*squeeze(1-VOC_perc_immun).*avg_trans_VOC_inf;
    ReffUK = (bet_UK./parameters.gam).*squeeze(1-UK_perc_immun).*avg_trans_UK_inf;
end

% Compute weighted force of infection, accounting for transmission blocking
% action of immunity
function [weighted_I_UK,weighted_I_VOC] = compute_weighted_force_of_infection(pop_out,p)
% Inputs: 
%  pop_out - proportion of popn in each state. Fourth dimension has slice per timestep
%  p - parameter structure
% Outputs: 
%  weighted_I_UK - Vector. Relative transmissibility of average infected
%                           with resident variants at each timestep.
%  weighted_I_VOC - Vector.Relative transmissibility of average infected
%                           with VOC MT each timestep.

    % Compute force of infection from unvaccinated and those suffering first
    % infection event
    I_UK_unvacc_first_infection = pop_out(3,1,1,:);
    I_VOC_unvacc_first_infection = pop_out(1,3,1,:);

    % Compute force of infection from unvaccinated infecteds who have had a prior infection
    I_UK_unvacc_prior_infection = pop_out(3,4,1,:)*p.propn_transmission_priorinf;
    I_VOC_unvacc_prior_infection = pop_out(4,3,1,:)*p.propn_transmission_priorinf;

    % Compute force of infection from infecteds (first infection event) who are vaccinated
    I_UK_vacc_first_infection = (pop_out(3,1,2,:)*p.propn_transmission_a) +... % AZ
                                 (pop_out(3,1,3,:)*p.propn_transmission_p) +... % Pfizer
                                 (pop_out(3,1,4,:)*p.propn_transmission_n); % new vaccine

    I_VOC_vacc_first_infection = (pop_out(1,3,2,:)*p.propn_transmission_a) +... % AZ
                                 (pop_out(1,3,3,:)*p.propn_transmission_p) +... % Pfizer
                                 (pop_out(1,3,4,:)*p.propn_transmission_n); % new vaccine

    % Compute force of infection from infecteds who are vaccinated AND had a prior infection
    prior_inf_AZ_min = min(p.propn_transmission_a,p.propn_transmission_priorinf);
    prior_inf_Pfizer_min = min(p.propn_transmission_p,p.propn_transmission_priorinf);
    prior_inf_newvacc_min = min(p.propn_transmission_n,p.propn_transmission_priorinf);

    I_UK_vacc_prior_infection = (pop_out(3,4,2,:)*prior_inf_AZ_min) +... % AZ
                                 (pop_out(3,4,3,:)*prior_inf_Pfizer_min) +... % Pfizer
                                 (pop_out(3,4,4,:)*prior_inf_newvacc_min); % new vaccine

    I_VOC_vacc_prior_infection = (pop_out(4,3,2,:)*prior_inf_AZ_min) +... % AZ
                                 (pop_out(4,3,3,:)*prior_inf_Pfizer_min) +... % Pfizer
                                 (pop_out(4,3,4,:)*prior_inf_newvacc_min); % new vaccine

    % Get cumulative force of infection.
    % Use squeeze to remove singleton dimensions and return a vector,
    % entry per timepoint
    weighted_I_UK = squeeze(I_UK_unvacc_first_infection + I_UK_unvacc_prior_infection +...
                    I_UK_vacc_first_infection + I_UK_vacc_prior_infection);
    weighted_I_VOC = squeeze(I_VOC_unvacc_first_infection + I_VOC_unvacc_prior_infection +...
                    I_VOC_vacc_first_infection + I_VOC_vacc_prior_infection);
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
%% plot_script_parsimonious_model.m:
% Produce figures of results generated using the parsimonious SARS-CoV-2
% transmission model and the stochastic importation model.
%--------------------------------------------------------------------------
clear
set(0,'defaultfigurecolor',[1 1 1])

%% check if the folder for saving figures exists, and if not then create it
if ~exist('saved_figs','dir')
    mkdir saved_figs
end

%% Add required directories to path
addpath('cbrewer/') % Colour maps

%% Set global flag variables
save_figs_flag = false;
glyph_flag = true;

%% Set global plot variables
% Set up plot markers
glyph_markertypes = {'s','+','o','v','*'};
glyph_markertypes_VOC_AtoC = {'s','+','o','none'};
glyph_markertypes_VOC_DtoE = {'+','v','*','none'};

% Colours to be used in line plots
line_colours_VOC_AtoC = [0    0.4470    0.7410;
                        0.8500    0.3250    0.0980;
                        0.9290    0.6940    0.1250;
                        0 0 0];

line_colours_VOC_DtoE = [0.8500    0.3250    0.0980;
                        0.4940    0.1840    0.5560;
                        0.4660    0.6740    0.1880;
                        0 0 0];                
                    
line_colours_new_vacc = [1 0 0;
                         1 0 0;
                         1 0 0;
                         1 0 0];

% Set line style type and linewidth for each trace/profile
line_type_VOC_AtoC = {'-','-','-'};
line_width_VOC_AtoC = [2 2 2];

line_type_VOC_DtoE = {'-','-','-'};
line_width_VOC_DtoE = [2 2 2];

line_type_VOC_AtoC_and_no_VOC = {'-','-','-','-'};
line_width_VOC_AtoC_and_no_VOC  = [2 2 2 2];

line_type_VOC_DtoE_and_no_VOC = {'-','-','-','-'};
line_width_VOC_DtoE_and_no_VOC  = [2 2 2 2];

line_type_new_vacc = {':',':',':',':'};
line_width_new_vacc = [2.5,2.5,2.5,2.5];

% Set y maximums
ymax_VOC_temporal = 3.8;
ymax_UK_temporal = 2.3;

ymax_VOC_Reff = 2.6;
ymax_UK_Reff = 2.6;

% Set cutoff for dates on xaxis in timeseries plots
xlim_end = datetime(2022,4,1);

%--------------------------------------------------------------------------
%% LOAD PARSIMONIOUS MODEL RESULTS DATA
%--------------------------------------------------------------------------
load('MAT_files/parsimonious_model_results.mat')

%--------------------------------------------------------------------------
%% VOC SCENARIO RUNS.
%% Temporal infectious prevalence with resident variants in absence of VOCs also displayed 
%% Figures 1(a), S5(a), S6
%--------------------------------------------------------------------------

% Get default parameter set. To be used in plotting function.
parameters = make_parameters();

%%% Figure including VOC MT, VOC E, VOC LT+E %%%
%%% Figures 1(a) & S5(a)

% Set up input data
UK_input_data = [I_UK_default_runs(:,1:3,1) I_UK_no_VOC];
VOC_input_data = [I_VOC_default_runs(:,1:3,1) I_UK_no_VOC];

% Set figure plot inputs
leg_labels = {'VOC MT','VOC E','VOC LT+E','Resident variants with no VOCs','January 2021 peak prevalence'};
leg_pos = [0.657 0.725 0.227 0.147];
save_filename_resident_variants = 'saved_figs/Temporal_resident_variants_and_no_VOC_trace';
save_filename_VOC =  'saved_figs/Temporal_VOCs_and_no_VOC_trace';
fig_fontsize = 22;

% Specify if traces from other VOC intro dates should be added
add_multiple_trace_flag = false;

% Specify if vaccination data should also be displayed
add_vaccination_data_flag = true;

% Specify if legend should be included
leg_flag = true;

% Generate figure
temporal_infectious_plots(outputs_default_runs,...
                          UK_input_data,...
                          VOC_input_data,...
                          I_UK_default_runs,...
                          I_VOC_default_runs,...
                          parameters,...
                          line_type_VOC_AtoC_and_no_VOC,...
                          line_colours_VOC_AtoC,...
                          line_width_VOC_AtoC_and_no_VOC,...
                          glyph_markertypes_VOC_AtoC,...
                          ymax_UK_temporal,...
                          ymax_VOC_temporal,...
                          xlim_end,...
                          leg_flag,...
                          leg_labels,...
                          leg_pos,...
                          save_figs_flag,...
                          save_filename_resident_variants,...
                          save_filename_VOC,...
                          fig_fontsize,...
                          [],...
                          add_multiple_trace_flag,...
                          vacc_coverage_data,...
                          add_vaccination_data_flag)
%%                      
%%% Figures show percentage of cases that are VOCs %%%
%%% Figures S7(a) and S7(b) %%%

% Set up input data
UK_input_data = I_UK_default_runs(:,:,1);
VOC_input_data = I_VOC_default_runs(:,:,1);

% Set figure plot inputs
leg_labels = {'VOC MT','VOC E','VOC LT+E','VOC Ev','VOC Ei'};
leg_pos = [0.657 0.725 0.227 0.147];
save_filename =  'saved_figs/Temporal_proportion_plot';
fig_fontsize = 22;

% Specify if vaccination data should also be displayed
add_vaccination_data_flag = true;

% Specify if legend should be included
leg_flag = true;

% Generate figure
temporal_proportion_plots(outputs_default_runs,...
                          I_UK_default_runs,...
                          I_VOC_default_runs,...
                          [line_type_VOC_AtoC_and_no_VOC; line_type_VOC_AtoC_and_no_VOC],...
                          [line_colours_VOC_AtoC(1:3,:); line_colours_VOC_DtoE(2:3,:)],...
                          [line_width_VOC_AtoC_and_no_VOC; line_width_VOC_AtoC_and_no_VOC],...
                          [glyph_markertypes_VOC_AtoC(1:3) glyph_markertypes_VOC_DtoE(2:3)],...
                          xlim_end,...
                          leg_flag,...
                          leg_labels,...
                          leg_pos,...
                          save_figs_flag,...
                          save_filename,...
                          fig_fontsize,...
                          vacc_coverage_data,...
                          add_vaccination_data_flag)
%%                      
%%% Figure including VOC Ev & VOC Ei %%%
%%% Figures S5(a) & S6(b)

% Set up input data
UK_input_data_batch_two = [I_UK_default_runs(:,[2 4:5],1) I_UK_no_VOC];
VOC_input_data_batch_two = [I_VOC_default_runs(:,[2 4:5],1) I_UK_no_VOC];

% Set figure plot inputs
leg_labels = {'VOC E','VOC Ev','VOC Ei','Resident variants with no VOCs','January 2021 peak prevalence'};
leg_pos = [0.657 0.725 0.227 0.147];
save_filename_resident_variants = 'saved_figs/Temporal_resident_variants_extra_and_no_VOC_trace';
save_filename_VOC =  'saved_figs/Temporal_VOCs_extra_and_no_VOC_trace';
fig_fontsize = 22;

% Specify if traces from other VOC intro dates should be added
add_multiple_trace_flag = false;

% Specify if vaccination data should be plotted
add_vaccination_data_flag = true;

% Specify if legend should be included
leg_flag = true;

% Generate figure
temporal_infectious_plots(outputs_default_runs,...
                          UK_input_data_batch_two,...
                          VOC_input_data_batch_two,...
                          I_UK_default_runs,...
                          I_VOC_default_runs,...
                          parameters,...
                          line_type_VOC_DtoE_and_no_VOC,...
                          line_colours_VOC_DtoE,...
                          line_width_VOC_DtoE_and_no_VOC,...
                          glyph_markertypes_VOC_DtoE,...
                          ymax_UK_temporal,...
                          ymax_VOC_temporal,...
                          xlim_end,...
                          leg_flag,...
                          leg_labels,...
                          leg_pos,...
                          save_figs_flag,...
                          save_filename_resident_variants,...
                          save_filename_VOC,...
                          fig_fontsize,...
                          [],...
                          add_multiple_trace_flag,...
                          vacc_coverage_data,...
                          add_vaccination_data_flag) 

                      
%--------------------------------------------------------------------------
%% VOC SCENARIO RUNS. R with immunity
%% Figures 1(c), S5(b)
%--------------------------------------------------------------------------                       
                      
%% Figure for VOCs MT, E, LT+E, Figure 1(c) %%%

% Set figure plot inputs
save_filename =  'saved_figs/R_with_immunity_default_VOCs';
fig_fontsize = 22;

% Generate figure
add_vaccination_data_flag = true;
temporal_R_with_immunity_plots(no_VOC_outputs,...
                          R0UK,...
                          [R0VOC{1:3}]',...
                          line_type_VOC_AtoC,...
                          line_colours_VOC_AtoC,...
                          line_width_VOC_AtoC,...
                          glyph_markertypes_VOC_AtoC,...
                          ymax_UK_Reff,...
                          save_figs_flag,...
                          save_filename,...
                          fig_fontsize,...
                          xlim_end,...
                          vacc_coverage_data,...
                          add_vaccination_data_flag) 

%% Figure for VOCs Ei, Ev, Figure S5(b) %%%
% Set figure plot inputs
save_filename =  'saved_figs/R_with_immunity_VOCs_extra';
fig_fontsize = 22;

% Generate figure
add_vaccination_data_flag = true;
temporal_R_with_immunity_plots(no_VOC_outputs,...
                          R0UK,...
                          [R0VOC{[2, 4:5]}]',...
                          line_type_VOC_DtoE,...
                          line_colours_VOC_DtoE,...
                          line_width_VOC_DtoE,...
                          glyph_markertypes_VOC_DtoE,...
                          ymax_UK_Reff,...
                          save_figs_flag,...
                          save_filename,...
                          fig_fontsize,...
                          xlim_end,...
                          vacc_coverage_data,...
                          add_vaccination_data_flag)

%--------------------------------------------------------------------------
%% SENSITIVITY HEATMAPS (FIGURES 1B & 1D) %%
%--------------------------------------------------------------------------

%% Change relative VOC_vs_UK transmissibility and AZ vaccine efficacy against VOC together

% Specify values to test
VOC_efficacy_scaling = 0.5:0.05:1;
VOC_vs_UK_varies = 0.5:0.1:1.5;

% Batch of single panel heatmaps. Sensitivity of VOC introduction date to other parameters

% Set plot properties
tick_label_fontsize = 22;
label_fontsize = 24;
glyph_pos = [11,11;
             6,6;
             4,6];
flip_yaxis_flag = true;

% Specify the save filename
save_filename_prefix = 'saved_figs/heatmap_immune_escape_vs_rel_trans_with_glyphs';

% Call function to construct plot
make_single_panel_heatmaps(VOC_vs_UK_varies,...
                            VOC_vs_UK_varies,...
                            VOC_efficacy_scaling(1:1:end),...
                            VOC_efficacy_scaling(1:1:end),...
                            VOC_vs_UK_varies,...
                            VOC_efficacy_scaling,...
                            epidemic_size_1(:,:,2),...
                            peak_height_1(:,:,2),...
                            peak_time_1(:,:,2),...
                            tick_label_fontsize,...
                            'Transmissibility of VOC vs resident variants',...
                            'Proportional efficacy against VOC',...
                            label_fontsize,'Variant of concern',...
                            base_outputs.date1,...
                            glyph_flag,...
                            glyph_pos,...
                            flip_yaxis_flag,...
                            save_figs_flag,...
                            save_filename_prefix)

                        
%%
%--------------------------------------------------------------------------
%% VOC E Gillespie part (FIGURE 3A) %%
%--------------------------------------------------------------------------
% Load data file
load('MAT_files/VOC_E_Gillespie')

% Set plot properties
label_fontsize = 26;
flip_yaxis_flag = true;
y_vals_is_date_type = false;

% plot over
% options: effective_imports_over, VOC_rel_trans_over,
% R_excl_immun_wildtype_over, relative_suscept_over
ordering = {'effective_imports_over','VOC_rel_trans_over','R_excl_immun_wildtype_over','relative_suscept_over'};
plot_over_x = 'effective_imports_over';
plot_over_y = 'VOC_rel_trans_over';

% Set axes labels
labels = {'Effective imports per day','Relative transmission','Resident variants R excluding immunity','relative immune escape'};
xlabel_string = labels(strcmp(ordering,plot_over_x));
ylabel_string = labels(strcmp(ordering,plot_over_y));

% set default values for other parameters
effective_imports_pos = find(effective_imports_over==0.2); % set effective imports per day
VOC_rel_trans_pos = find(VOC_rel_trans_over==1);
R_excl_immun_wildtype_pos = find(R_excl_immun_wildtype_over==3);
relative_suscept_pos = 1;
Index = {effective_imports_pos,VOC_rel_trans_pos,R_excl_immun_wildtype_pos,relative_suscept_pos};

% find matrix for plot_over_x and y
Index{strcmp(ordering,plot_over_x)} = ':';
Index{strcmp(ordering,plot_over_y)} = ':';
matrix_to_plot = squeeze(epidemic_prob(Index{:}))';
if find(strcmp(ordering,plot_over_x))>find(strcmp(ordering,plot_over_y))
    matrix_to_plot = matrix_to_plot';
end

% Set up colourmap
CT_map = cbrewer('seq','Oranges', 128);

% Set up colourbar range
cbar_range = [0 1];

% plot epidemic probability
title_string = 'Epidemic probability';
cbar_type = 4;
save_filename = 'saved_figs/epidemic_probability';
make_heatmap_plot_new_vacc(matrix_to_plot,...
    eval(plot_over_x),...
    eval(plot_over_y),...
    label_fontsize,...
    xlabel_string,...
    ylabel_string,...
    title_string,...
    flip_yaxis_flag,...
    y_vals_is_date_type,...
    cbar_type,...
    cbar_range,...
    CT_map,...
    save_filename)

%%
%--------------------------------------------------------------------------
%% NEW VACCINE VOC SUMMARY STAT HEATMAPS - RUNS E (FIGURE 3B-D)      %%
%% Previously vaccinated prioritised to receive VOC-targeted vaccine %%
%--------------------------------------------------------------------------

%  Epidemic size VOC, peak height VOC & peak time VOC

% Load data file
load('MAT_files/VOC_E_runs_vacc_priority_5')

% Set plot inputs
z_value = 1;
label_fontsize = 26;
flip_yaxis_flag = true;
y_vals_is_date_type = true;

% Set up axes labels
plot_over_x_name = 'Effective imports per day';
plot_over_z_name = 'Date VOC targeted vaccine introduced';

% Set colourbar ranges
% Row 1: Minimum; Row 2: Maximum
% Column for each summary statistics: Size, peak height, peak time
% Slice for resident variants, VOC, both
cbar_ranges = zeros(2,3,3);
cbar_ranges(:,:,1) = [40 0.5 datenum(2021,8,1);
                      50 1 datenum(2022,1,1)];

cbar_ranges(:,:,2) = [0 0 datenum(2021,6,1);
                      50 2.0 datenum(2022,5,1)];  
                  
cbar_ranges(:,:,3) = [80 0.95 datenum(2021,8,1);
                      85.5 2.4 datenum(2022,1,1)];

% Set up colourmap
CT_map = cbrewer('seq', 'Greens', 128);

% Set save filenames
save_filename_prefix = 'saved_figs/';
save_filename_variant_type = {'resident_variants_run_E_vacc_priority_5',...
                              'VOC_run_E_vacc_priority_5',...
                              'both_run_E_vacc_priority_5'};

% Produce plots
for itr = 2:2
    % plot epidemic size
    title_string = 'Attack rate';
    if itr==1
        title_string = [title_string,' resident variants (% of population)'];
    elseif itr==2
        title_string = [title_string,' VOC (% of population)'];
    else
        title_string = [title_string,' both (% of population)'];
    end
    cbar_type = 3;
    save_filename = [save_filename_prefix 'epidemic_size_' save_filename_variant_type{itr} '_26May2021'];
    make_heatmap_plot_new_vacc(squeeze(median(epidemic_size_gillespie(6:end,z_value,:,:,itr),4)*100)',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        y_vals_is_date_type,...
        cbar_type,...
        cbar_ranges(:,1,itr),...
        CT_map,...
        save_filename)
    
    % plot peak height
    title_string = 'Peak in infectious prevalence for';
    if itr==1
        title_string = [title_string,' resident variants (% of population)'];
    elseif itr==2
        title_string = [title_string,' VOC (% of population)'];
    else
        title_string = [title_string,' both (% of population)'];
    end
    cbar_type = 3;
    save_filename = [save_filename_prefix 'peak_height_' save_filename_variant_type{itr} '_26May2021'];
    make_heatmap_plot_new_vacc(squeeze(median(peak_height_gillespie(6:end,z_value,:,:,itr),4)*100)',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        y_vals_is_date_type,...
        cbar_type,...
        cbar_ranges(:,2,itr),...
        CT_map,...
        save_filename)
    
    % plot peak time
    title_string = 'Time of infectious prevalence peak for';
    if itr==1
        title_string = [title_string,' resident variants'];
    elseif itr==2
        title_string = [title_string,' VOC'];
    else
        title_string = [title_string,' both'];
    end
    cbar_type = 2;
    save_filename = [save_filename_prefix 'peak_time_' save_filename_variant_type{itr} '_26May2021'];
    make_heatmap_plot_new_vacc(squeeze(median(datenum(peak_time_gillespie(6:end,z_value,:,:,itr)),4))',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        y_vals_is_date_type,...
        cbar_type,...
        cbar_ranges(:,3,itr),...
        CT_map,...
        save_filename)
end                        

%%
%--------------------------------------------------------------------------
%% NEW VACCINE VOC SUMMARY STAT HEATMAPS - RUNS E (FIGURE S13)      %%
%% Previously vaccinated prioritised to receive VOC-targeted vaccine %%
%--------------------------------------------------------------------------

%  Epidemic size VOC, peak height VOC & peak time VOC

% Load data file
load('MAT_files/VOC_E_runs_vacc_priority_4')

% Set plot inputs
z_value = 1;
label_fontsize = 26;
flip_yaxis_flag = true;
y_vals_is_date_type = true;

% Set up axes labels
plot_over_x_name = 'Effective imports per day';
plot_over_z_name = 'Date VOC targeted vaccine introduced';

% Set colourbar ranges
% Row 1: Minimum; Row 2: Maximum
% Column for each summary statistics: Size, peak height, peak time
% Slice for resident variants, VOC, both
cbar_ranges = zeros(2,3,3);
cbar_ranges(:,:,1) = [40 0.5 datenum(2021,8,1);
                      50 1 datenum(2022,1,1)];

cbar_ranges(:,:,2) = [0 0 datenum(2021,6,1);
                      50 2.0 datenum(2022,5,1)];  
                  
cbar_ranges(:,:,3) = [80 0.95 datenum(2021,8,1);
                      85.5 2.4 datenum(2022,1,1)];

% Set up colourmap
CT_map = cbrewer('seq', 'Greens', 128);

% Set save filenames
save_filename_prefix = 'saved_figs/';
save_filename_variant_type = {'resident_variants_run_E_vacc_priority_4',...
                              'VOC_run_E_vacc_priority_4',...
                              'both_run_E_vacc_priority_4'};

% Produce plots
for itr = 2:2
    % plot epidemic size
    title_string = 'Attack rate';
    if itr==1
        title_string = [title_string,' resident variants (% of population)'];
    elseif itr==2
        title_string = [title_string,' VOC (% of population)'];
    else
        title_string = [title_string,' both (% of population)'];
    end
    cbar_type = 3;
    save_filename = [save_filename_prefix 'epidemic_size_' save_filename_variant_type{itr} '_26May2021'];
    make_heatmap_plot_new_vacc(squeeze(median(epidemic_size_gillespie(6:end,z_value,:,:,itr),4)*100)',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        y_vals_is_date_type,...
        cbar_type,...
        cbar_ranges(:,1,itr),...
        CT_map,...
        save_filename)
    
    % plot peak height
    title_string = 'Peak in infectious prevalence for';
    if itr==1
        title_string = [title_string,' resident variants (% of population)'];
    elseif itr==2
        title_string = [title_string,' VOC (% of population)'];
    else
        title_string = [title_string,' both (% of population)'];
    end
    cbar_type = 3;
    save_filename = [save_filename_prefix 'peak_height_' save_filename_variant_type{itr} '_26May2021'];
    make_heatmap_plot_new_vacc(squeeze(median(peak_height_gillespie(6:end,z_value,:,:,itr),4)*100)',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        y_vals_is_date_type,...
        cbar_type,...
        cbar_ranges(:,2,itr),...
        CT_map,...
        save_filename)
    
    % plot peak time
    title_string = 'Time of infectious prevalence peak for';
    if itr==1
        title_string = [title_string,' resident variants'];
    elseif itr==2
        title_string = [title_string,' VOC'];
    else
        title_string = [title_string,' both'];
    end
    cbar_type = 2;
    save_filename = [save_filename_prefix 'peak_time_' save_filename_variant_type{itr} '_26May2021'];
    make_heatmap_plot_new_vacc(squeeze(median(datenum(peak_time_gillespie(6:end,z_value,:,:,itr)),4))',...
        plot_over_x(6:end),...
        plot_over_z,...
        label_fontsize,...
        plot_over_x_name,...
        plot_over_z_name,...
        title_string,...
        flip_yaxis_flag,...
        y_vals_is_date_type,...
        cbar_type,...
        cbar_ranges(:,3,itr),...
        CT_map,...
        save_filename)
end           

%--------------------------------------------------------------------------
%% SENSITIVITY TO INTRODUCTION DATE CLOUD PLOTS %%
%% Figure S8
%--------------------------------------------------------------------------
 
%%% Cloud plots comparing infections, peak time and peak height %%%

% Set figure plot inputs
leg_labels = {'VOC MT','VOC E','VOC LT+E'};
leg_pos = [0.67 0.75 0.15 0.147];
fig_fontsize = 24;

% Set save filenames
save_filename_resident_variants = 'saved_figs/default_VOC_runs/cloud_plots_resident_variants';
save_filename_VOC =  'saved_figs/default_VOC_runs/cloud_plots_VOC';

% Call plot function
generate_summary_measure_cloud_plots(epidemic_size_default_runs(:,1:3,:),...
                                    peak_height_default_runs(:,1:3,:),...
                                    peak_time_default_runs,...
                                    line_colours_VOC_AtoC,...
                                    glyph_markertypes_VOC_AtoC,...
                                    leg_labels,...
                                    leg_pos,...
                                    fig_fontsize,...
                                    save_filename_resident_variants,...
                                    save_filename_VOC,...
                                    save_figs_flag)                              
                                
%--------------------------------------------------------------------------
%% Heatmaps of effective R as transmission and immune escape varies
%% Figure S9
%--------------------------------------------------------------------------

% Construct single panel heatmaps for relative effective R > 1 & R_VOC > 1
% Made as snapshots in time for different calendar dates.
% Date for seeding initial VOC infecteds & time snapshots to compute
% relative R
VOC_imp_date_varies = [datenum(2021,5,17) datenum(2021,8,1) datenum(2021,11,1)];
t_snapshots = [datenum(2021,5,17) datenum(2021,8,1) datenum(2021,11,1)];

% Transmissibility varied
VOC_vs_UK_varies = 0.5:0.1:1.5;

% Vaccine efficacy varied
VOC_efficacy_scaling = 0.5:0.05:1;

% Store outputs & labels in cells for use in plotting
save_filename_suffix = {'rel_eff_R','eff_R_VOC','eff_R_resident_variants','rel_eff_R_threshold'};

% Set plot properties
tick_label_fontsize = 22;
label_fontsize = 22;

% Get relevant data from cell store variable
R_eff_heatmap = effective_R_data{4};

% Populate each panel
for panel_idx = 1:numel(t_snapshots)

    % Initialise the figure
    position = [10, 10, 1.2*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(panel_idx);
    clf;
    set(fig,'Color', [1 1 1])
    hold on

    % Produce the heatmap
    imagesc(VOC_vs_UK_varies,VOC_efficacy_scaling,R_eff_heatmap(:,:,panel_idx));
    set(gca,'YDir','normal'); % Reinvert y-axis
    
    % Axis tick labels
    xticks(VOC_vs_UK_varies); xticklabels(VOC_vs_UK_varies);
    yticks(VOC_efficacy_scaling(1:1:end));
    yticklabels(VOC_efficacy_scaling(1:1:end));
    
    % Set title
    title(['Date: ', datestr(t_snapshots(panel_idx))]);
    
    % Add axis labels if required
    xlabel('Transmissibility of VOC vs resident variants','FontSize',label_fontsize)
    ylabel('Proportional efficacy against VOC','FontSize',label_fontsize)
    
    % Plot properties
    set(gca,'FontSize',tick_label_fontsize);
    
    % Set ticklabel properties
    ytickformat('%.2f')
    xtickformat('%.1f')
    
    % Add letters
    glyph_pos = [11,11;
        6,6;
        4,6];
    glyph_size = 18;
    glyph_markertypes = {'MT','E','LT+E'};
    glyph_colour = [1 0 0];
    if glyph_flag == true
        add_glyph_fn(glyph_pos,glyph_markertypes,glyph_colour,glyph_size,...
            VOC_vs_UK_varies,... % xticks values
            VOC_efficacy_scaling(1:1:end)) %yticks values
    end
    
    % Set colourmap
    colormap([1 1 1;
            0 0 0])
        
    % Set axes limits
    xlim([0.48 1.52])
    ylim([0.48 1.02])

    % Set figure properties
    % axis image
    box on
        
    % save figure
    if save_figs_flag == true
        save_filename = ['saved_figs/relative_R/heatmap_', save_filename_suffix{4},'_panel_idx_' num2str(panel_idx)];
        export_fig(save_filename,'-pdf','-painters','-r1200')
    end
end

%--------------------------------------------------------------------------
%% Heatmaps of effective R as transmission and immune escape varies
%% Extra figures
%--------------------------------------------------------------------------
% Construct array of heatmaps for 
% (i) relative effective R
% (ii) R_VOC
% (iii) R_resident_variants
% (iv) relative effective R > 1 & R_VOC > 1

% Store labels in cells for use in plotting
colourbar_label = {'Relative effective R (R_{eff}^{VOC}/R_{eff}^{res})','VOC effective R, R_{eff}^{VOC}','Resident variants effective R, R_{eff}^{res}','Relative effective R > 1'};
save_filename_suffix = {'rel_eff_R','eff_R_VOC','eff_R_resident_variants','rel_eff_R_threshold'};

% Set plot properties
tick_label_fontsize = 22;
label_fontsize = 22;

% Constuct the heatmaps
for fig_itr = 1:4
    position = [10, 10, 3.5*550, 1.2*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(fig_itr);
    clf;
    set(fig,'Color', [1 1 1])
    hold on

    % Get relevant data from cell store variable
    R_eff_heatmap = effective_R_data{fig_itr};

    % Populate each panel
    for panel_idx = 1:numel(t_snapshots)
        subplot(1,3,panel_idx)
        
        % Produce the heatmap
        CT_purple = cbrewer('seq', 'Purples', 128);
        imagesc(VOC_vs_UK_varies,VOC_efficacy_scaling,R_eff_heatmap(:,:,panel_idx));
        set(gca,'YDir','normal'); % Reinvert y-axis

        % Add colourbar
        if fig_itr ~= 4
            c = colorbar;
            c.FontSize = 20;
            if mod(panel_idx,3) == 0
                hL = ylabel(c,colourbar_label{fig_itr});
                %set(hL,'Rotation',270);
            end
        end
        
        % Set colourbar range
        c_max = max(R_eff_heatmap(:));
        caxis([0,c_max]);

        % Axis tick labels
        xticks(VOC_vs_UK_varies); xticklabels(VOC_vs_UK_varies);
        yticks(VOC_efficacy_scaling(1:1:end)); 
        yticklabels(VOC_efficacy_scaling(1:1:end));
        %yticks(efficacy_varies(end:-2:1)); 
        %yticklabels(efficacy_varies(1:2:end));

        % Set title 
        title(['Date: ', datestr(t_snapshots(panel_idx))]);
        
        % Add axis labels if required
        if fig_itr == 4
            xlabel('Transmissibility of VOC vs resident variants','FontSize',label_fontsize)
            ylabel('Proportional efficacy against VOC','FontSize',label_fontsize)
        else
            if panel_idx == 2
                xlabel('Transmissibility of VOC vs resident variants','FontSize',label_fontsize)
            end
            if mod(panel_idx,3) == 1
                ylabel('Proportional efficacy against VOC','FontSize',label_fontsize)
                %ylabel('Susceptibility to VOC post-AZ vaccine (e_{aVOC})','FontSize',label_fontsize);
            end
        end
        hold on
        
        % Plot properties
        set(gca,'FontSize',tick_label_fontsize);
        
        % Set ticklabel properties
        ytickformat('%.2f')
        xtickformat('%.1f')
        
        % Add letters
        glyph_pos = [11,11;
                     6,6;
                     4,6];
        glyph_size = 18;
        glyph_markertypes_vec = {'MT','E','LT+E'};  
        glyph_colour = [1 0 0];
        if glyph_flag == true
            add_glyph_fn(glyph_pos,glyph_markertypes_vec,glyph_colour,glyph_size,...
                VOC_vs_UK_varies,... % xticks values
                VOC_efficacy_scaling(1:1:end)) %yticks values
        end
        
        % Set colourmap
        if fig_itr ~= 4
            colormap(CT_purple)
        elseif fig_itr == 4
            colormap([1 1 1;
                      0 0 0])
        end  
    end
    
    % save figure
    if save_figs_flag == true
        save_filename = ['saved_figs/relative_R/heatmap_', save_filename_suffix{fig_itr}];
        export_fig(save_filename,'-pdf','-painters','-r1200')
    end
end

%--------------------------------------------------------------------------
%% TRANSMISSION BLOCKING SENSITIVITY RUNS  %%
%% Figures S10 & S11
%--------------------------------------------------------------------------

% Get default parameter set. To be used in plotting function.
parameters = make_parameters();

% Options for varying the transmission blocking assumptions
transmission_blocking_val = [0 0.25 0.5];
n_transmission_blocking_vals = numel(transmission_blocking_val);

% Set up VOC parameters
VOC_vs_UK_varies = [1.5 1   0.8]; % Relative transmissibility of VOC
efficacy_varies =  [1   0.75 0.75]; % Vaccine efficacy
s_varies =         [1   0.75 0.75]; % Natural immunity efficacy
n_VOCs = numel(VOC_vs_UK_varies);

%% Generate temporal plots
%% Figure S10

% Set figure plot inputs
leg_labels = {'VOC MT','VOC E','VOC LT+E','Resident variants with no VOCs','January 2021 peak prevalence'};
leg_pos = [0.657 0.725 0.227 0.147];
fig_fontsize = 23;
add_vaccination_data_flag = true;
add_multiple_trace_flag = false;

% Generate figure
for jj = 2:n_transmission_blocking_vals
    
    if jj == 2
        leg_flag = true;
    else
        leg_flag = false;
    end
    
    % Set up input data
    UK_input_data_batch_trans_block = [I_UK_trans_block(:,1:3,jj) I_UK_no_VOC_trans_block(:,jj)];
    VOC_input_data_batch_trans_block = [I_VOC_trans_block(:,1:3,jj) I_UK_no_VOC_trans_block(:,jj)];

    % Set function inputs
    save_filename_resident_variants = ['saved_figs/transmission_blocking/Temporal_resident_variants_transmission_blocking_', num2str(transmission_blocking_val(jj)*100),'percent'];
    save_filename_VOC =  ['saved_figs/transmission_blocking/Temporal_VOC_transmission_blocking_', num2str(transmission_blocking_val(jj)*100),'percent'];
    title_input = ['Transmission blocking: ',num2str(transmission_blocking_val(jj)*100),'%'];
    
    % Call plot function
    temporal_infectious_plots(outputs_trans_block,...
                              UK_input_data_batch_trans_block,...
                              VOC_input_data_batch_trans_block,...
                              I_UK_trans_block,...
                              I_VOC_trans_block,...
                              parameters,...
                              line_type_VOC_AtoC_and_no_VOC,...
                              line_colours_VOC_AtoC,...
                              line_width_VOC_AtoC_and_no_VOC,...
                              glyph_markertypes_VOC_AtoC,...
                              ymax_UK_temporal,...
                              ymax_VOC_temporal,...
                              xlim_end,...
                              leg_flag,...
                              leg_labels,...
                              leg_pos,...
                              save_figs_flag,...
                              save_filename_resident_variants,...
                              save_filename_VOC,...
                              fig_fontsize,...
                              title_input,...
                              add_multiple_trace_flag,...
                              vacc_coverage_data,...
                              add_vaccination_data_flag)   
end

%% Generate side-by-side effective R plots for each variant
%% Figure S11

% Array for colours to be used in plot
line_colours = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560];

% Set line style type and linewidth for each trace/profile
line_type = {'-','-','-'};
line_width = [2 2 2];

% Get timestep index to plot from for VOC
VOC_introduced_timestep_idx = (parameters.VOC_imp_date - parameters.date1) + 1;

% Set x tick values to be used in upcoming plots
xtick_vals = datetime(2021,6,1) + calmonths(0:2:11);
    
% Initialise figure
position = [10, 10, 3.2*550, 2.1*450];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])
hold on
    
% Populate each panel
% Column per transmission blocking value.
% Row 1 for resident variants. Row 2 for VOC.
for jj = 1:n_transmission_blocking_vals

    for VOC_itr = 1:n_VOCs
        subplot(2,n_transmission_blocking_vals,jj)
        hold on
        % Get effective R data from ReffVOC_trans_block
        % First dimension for transmission blocking value.
        % Second dimension for each VOC.
        % Third dimension for time.
        UK_data_to_plot = squeeze(ReffUK_trans_block(jj,VOC_itr,:));

        % Plot outputs
        plot(outputs_trans_block.dates,...
                        UK_data_to_plot,...
                        line_type{VOC_itr},...
                        'Color',line_colours(VOC_itr,:),...
                        'LineWidth',line_width(VOC_itr),...
                        'Marker',glyph_markertypes{VOC_itr},...
                        'MarkerSize',12,...
                        'MarkerIndices',1:20:length(UK_data_to_plot));
        
        subplot(2,n_transmission_blocking_vals,n_transmission_blocking_vals+jj)
        hold on
        % Get effective R data from ReffUK_trans_block
        % First dimension for transmission blocking value.
        % Second dimension for each VOC.
        % Third dimension for time.
        VOC_data_to_plot = squeeze(ReffVOC_trans_block(jj,VOC_itr,VOC_introduced_timestep_idx:end));

        % Plot outputs
        plot(outputs_trans_block.dates(VOC_introduced_timestep_idx:end),...
                        VOC_data_to_plot,...
                        line_type{VOC_itr},...
                        'Color',line_colours(VOC_itr,:),...
                        'LineWidth',line_width(VOC_itr),...
                        'Marker',glyph_markertypes{VOC_itr},...
                        'MarkerSize',12,...
                        'MarkerIndices',1:20:length(VOC_data_to_plot));            
    end
    
    % Subplot figure properties
    subplot(2,n_transmission_blocking_vals,jj)
    % x-axis properties
    xlim([outputs_trans_block.dates(1),outputs_trans_block.dates(end)])
    xticks(xtick_vals)
    
    % y-axis properties
    ylim([0 2.5])
    if jj == 1
        ylabel('Effective R for resident variants, R_{eff}^{res}');
    end
    
    % Add title
    title(['Transmission blocking: ',num2str(transmission_blocking_val(jj)*100),'%'],...
               'FontWeight','bold',...
               'FontSize',22)
           
    % Add legend
    if jj == n_transmission_blocking_vals
        legend({'VOC MT','VOC E','VOC LT+E'},...
                'FontSize',22)
    end
    
    % Set plot properties
    set(gca,'Fontsize',22)
    set(gca,'LineWidth',1)
    box on
    
    subplot(2,n_transmission_blocking_vals,n_transmission_blocking_vals+jj)
    % x-axis properties
    xlabel('Time');
    xlim([outputs_trans_block.dates(1),outputs_trans_block.dates(end)])
    xticks(xtick_vals)
    
    % y-axis properties
    ylim([0 2.5])
    if jj == 1
        ylabel('Effective R for VOC, R_{eff}^{VOC}')
    end
    
    % Set plot properties
    set(gca,'Fontsize',22)
    set(gca,'LineWidth',1)
    box on
end

% Save figure
if save_figs_flag == true
    export_fig('saved_figs/transmission_blocking/R_eff_comparison','-pdf','-painters','-r1200')
end
                      
%--------------------------------------------------------------------------
%% Supporting functions to process outputs %%
%--------------------------------------------------------------------------

% Plot VOC scenario temporal plots
function temporal_infectious_plots(outputs,...
                                    I_UK,...
                                    I_VOC,...
                                    I_UK_alt_intro_dates,...
                                    I_VOC_alt_intro_dates,...
                                    parameters,...
                                    line_type,...
                                    line_colours,...
                                    line_width,...
                                    glyph_markertypes,...
                                    ymax_UK_temporal,...
                                    ymax_VOC_temporal,...
                                    xlim_end,...
                                    leg_flag,...
                                    leg_labels,...
                                    leg_pos,...
                                    save_figs_flag,...
                                    save_filename_resident_variants,...
                                    save_filename_VOC,...
                                    fig_fontsize,...
                                    title_input,...
                                    add_multiple_trace_flag,...
                                    vacc_data,...
                                    add_vaccination_data_flag)

    % Get timestep index to plot from for VOC
    VOC_introduced_timestep_idx = (parameters.VOC_imp_date - parameters.date1) + 1;

    % Set x tick values to be used in upcoming plots
    xtick_vals = datetime(2021,6,1) + calmonths(0:2:11);

    % Figure for VOCs
    for fig_itr = 1:2
        position = [10, 10, 2.8*550, 1.5*450];
        set(0, 'DefaultFigurePosition', position);
        fig = figure(fig_itr);
        clf;
        set(fig,'Color', [1 1 1])
        hold on

        % Add to figure. Output depends on index of loop.
        if fig_itr == 1
            % Resident variants
            y = I_UK(:,:,1)*100;
            y_incl_other_intro_dates = I_UK_alt_intro_dates*100;
            x_start_idx = 1;
        else
            % VOC
            x_start_idx = VOC_introduced_timestep_idx;
            y = I_VOC(x_start_idx:end,:,1)*100;
            y_incl_other_intro_dates = I_VOC_alt_intro_dates(x_start_idx:end,:,:)*100;
        end
        x = outputs.dates(x_start_idx:end);

        % Specify plotting against left y-axis
        yyaxis left
        ax_left = gca;

        % Add profiles to figure for each VOC
        for VOC_itr = 1:size(y,2)
            
            % Set colour for line
            if fig_itr == 1
                % Resident variants
                line_colour_vec = [0 0 0];
            else
                % VOC
                line_colour_vec = line_colours(VOC_itr,:);
            end

            % Add line profile for VOC
            plot(x,...
                    y(:,VOC_itr),...        
                    line_type{VOC_itr},...
                    'Color',line_colour_vec,...
                    'LineWidth',line_width(VOC_itr),...
                    'Marker',glyph_markertypes{VOC_itr},...
                    'MarkerSize',12,...
                    'MarkerIndices',1:20:length(y(:,VOC_itr)))

            % Add intervals for a range of introduction dates, if applicable
            if add_multiple_trace_flag == true
                % Data from other VOC intro dates stored in slices of y_incl_other_intro_dates
                for trace_itr = 2:size(y_incl_other_intro_dates,3)
                   h =  plot(x,...
                        y_incl_other_intro_dates(:,VOC_itr,trace_itr),...
                        line_type{VOC_itr},...
                        'Color',rgba2rgb([1 1 1],[line_colours(VOC_itr,:) 0.5]),...
                        'LineWidth',line_width(VOC_itr),...
                        'DisplayName','');

                    % Do not include in legend
                   set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
            end
        end
        plot(outputs.dates,ones(size(outputs.dates))*2,'k--','LineWidth',1)

        % y-axis properties
        if fig_itr == 1
            ylim([0 ymax_UK_temporal])
        else
            ylim([0 ymax_VOC_temporal])
        end

        % Set axis colour 
       ax_left.YAxis(1).Color = [0 0 0];

       % Panel specific
        if fig_itr == 1
            ylabel('Resident variants infectious prevalence');

            % Add title (if applicable)
            if ~isempty(title_input)
                title(title_input,...
                    'FontWeight','bold',...
                    'FontSize',22,...
                    'Units', 'normalized', 'Position', [0.5, 1.01, 0])
            end
        else
            ylabel('VOC infectious prevalence')
        end

        % Set y-axis tick label style
        ytickformat('%.1f%%')

        % x-axis properties
        xlabel('Time');
        xlim([outputs.dates(1),xlim_end])
        xticks(xtick_vals)

        % Add vaccination data if required
        if add_vaccination_data_flag == true

            if fig_itr == 1
                resident_variant_plot_flag = true;
            else
                resident_variant_plot_flag = false;
            end

            % Add vaccination onto the same panel
            add_vaccination_data_same_panel(outputs.dates,...
                vacc_data,...
                fig_fontsize,...
                resident_variant_plot_flag)
        end

        % Add legend
        if (leg_flag == true) %&& (fig_itr == 1)
            leg = legend(leg_labels,...
                            'Position',leg_pos);

            % Append vaccination bar data descriptor to legend
            bar_colour = rgba2rgb([1 1 1],[0 0 0.8 0.3]);
            bar(NaT,NaN,...
                            'FaceColor',bar_colour,...
                            'EdgeColor',bar_colour,...
                            'LineWidth',0.01,...
                            'DisplayName','Vaccination uptake (%)');

            % Set up indexing to move final label to be first, and shift all other
            % labels down by one
            n_labels = numel(leg.PlotChildren);
            label_reorder_idx = zeros(n_labels,1);
            label_reorder_idx(1) = n_labels;
            label_reorder_idx(2:end) = 1:(n_labels-1);

            % Reorder the legend labels
            reorder_leg = leg.PlotChildren(label_reorder_idx);
            legend(reorder_leg,'Position',leg_pos);
        end

        % Set plot properties
        set(ax_left,'Fontsize',fig_fontsize)
        set(ax_left,'LineWidth',1)
        box on

        % Save figures
        if save_figs_flag == true
            if fig_itr == 1
                export_fig(save_filename_resident_variants,'-pdf','-painters','-r1200')
            else
                export_fig(save_filename_VOC,'-pdf','-painters','-r1200')
            end
        end
    end 
end

% Add vaccination data to infectious prevalence timeseries plots
% Also add vertical lines denoting roadmap relaxation steps
function add_vaccination_data_same_panel(x_axis_dates,...
                                        vacc_data,...
                                        fig_fontsize,...
                                        resident_variant_plot_flag)
                       
    % Add data to right-hand y-axis 
    yyaxis right
    ax2 = gca;
    bar_colour = [0.0 0.0 0.8];
    alpha_val = 0.1;
    bar1 = bar(x_axis_dates,vacc_data,...
                'FaceColor',bar_colour,...
                'EdgeColor',bar_colour,...,
                'LineWidth',0.01);           
    bar1.FaceAlpha = alpha_val;
    bar1.EdgeAlpha = alpha_val;
    %set(ax2,'XAxisLocation','top','YAxisLocation','right','ydir','reverse');

    % Set y-axis labels
    yticks([0 20 40 60 80 100])
    ylim([0 105])

    % Set y-axis tick label style and label
    ytickformat('%g%%')
    ylabel('Vaccinated')

    % Set axis colour
    bar_colour_alpha = rgba2rgb([1 1 1],[bar_colour 0.5]);
    ax2.YAxis(2).Color = bar_colour_alpha;

    % Set height to add labels for each step of relaxation
    if resident_variant_plot_flag == true
        % Plots for resident variants
        y_text_pos = 90;
    else
        % Plots for VOC
        y_text_pos = 99;
    end

    % Add vertical lines for relaxation steps
    % Add STEP 3 date line
    text(datenum(2021,5,17)+3-datenum(x_axis_dates(1)),y_text_pos,{'Step 3';'R_{excl} = 2.41'},'Rotation',0,'HorizontalAlignment','Left','VerticalAlignment','Top','FontSize',fig_fontsize);

    % Add STEP 4 date line
    plot(datetime(2021,6,21)+[0 0],[0 105],'-',...
            'Color',[0.5 0.5 0.5],...
            'LineWidth',1);
    text(datetime(2021,6,21)+3,y_text_pos,{'Step 4';'R_{excl} = 3.51'},'Rotation',0,'HorizontalAlignment','Left','VerticalAlignment','Top','FontSize',fig_fontsize);
end


% Function to plot R with immunity %%
function  temporal_R_with_immunity_plots(outputs_default_runs,...
                          ReffUK_default_run,...
                          ReffVOC_default_runs,...
                          line_type,...
                          line_colours,...
                          line_width,...
                          glyph_markertypes,...
                          ymax_UK_Reff,...
                          save_figs_flag,...
                          save_filename,...
                          fig_fontsize,...
                          xlim_end,...
                          vacc_data,...
                          add_vaccination_data_flag)                     
    % Set x tick values to be used in upcoming plots
    xtick_vals = datetime(2021,6,1) + calmonths(0:2:11);

    % Initialise figures
    position = [10, 10, 2.8*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure;
    clf;
    set(fig,'Color', [1 1 1])
    hold on

    for VOC_itr = 1:size(ReffVOC_default_runs,1)
        % Get effective R data from ReffVOC_trans_block
        % First dimension for transmission blocking value.
        % Second dimension for each VOC.
        % Third dimension for time.
        UK_data_to_plot = squeeze(ReffVOC_default_runs(VOC_itr,:));

        % Plot outputs
        plot(outputs_default_runs.dates,...
            UK_data_to_plot,...
            line_type{VOC_itr},...
            'Color',line_colours(VOC_itr,:),...
            'LineWidth',line_width(VOC_itr),...
            'Marker',glyph_markertypes{VOC_itr},...
            'MarkerSize',12,...
            'MarkerIndices',1:20:length(UK_data_to_plot));
    end

    % plot the resident variants
    UK_data_to_plot = ReffUK_default_run;

    % Plot outputs
    plot(outputs_default_runs.dates,UK_data_to_plot,'k-','LineWidth',line_width(1));

    % Add labels etc to each figure
    % x-axis properties
    xlabel('VOC introduction date');
    xlim([outputs_default_runs.dates(1),xlim_end])
    xticks(xtick_vals)

    % y-axis properties
    ylabel('R with immunity');
    ylim([0.6 ymax_UK_Reff])

    % Add vaccination data if required
    if add_vaccination_data_flag == true
        resident_variants_plot_flag = false;    
        add_vaccination_data_same_panel(outputs_default_runs.dates,...
                                        vacc_data,...
                                        fig_fontsize,...
                                        resident_variants_plot_flag)
    end

    % Set plot properties
    set(gca,'Fontsize',fig_fontsize)
    set(gca,'LineWidth',1)
    box on

    if ispc
        g = gca;
        g.FontSize = g.FontSize*72/96;
    end

    % Save figures
    if save_figs_flag == true
        export_fig(save_filename,'-pdf','-painters','-r1200')
    end
end


% Batch of three heatmaps. Each saved individually
% Epidemics size; peak size; peak timing.
function make_single_panel_heatmaps(xticks_vals,xticks_labels,yticks_vals,yticks_labels,...
                                    x_data,y_data,...
                                    epidemic_size,peak_height,peak_time,...
                                    tick_label_fontsize,xlabel_text,ylabel_text,...
                                    label_fontsize,variant,date1,...
                                    glyph_flag,glyph_pos,flip_yaxis_flag,...
                                    save_figs_flag,save_filename_prefix)

    if glyph_flag == true
        % Set markersize, markertypes and colour values
        glyph_size = 18;
        % glyph_markertypes = {'A','B','C','D'};
        glyph_markertypes = {'MT','E','LT+E'};  
        glyph_colour = [1 0 0];
    end

    % Set colour maps
    CT_purple = cbrewer('seq', 'Purples', 128);

    % Final size
    f = figure;
    f.Position = [ 3   40   560*1.2   420*1.5];
    t = tiledlayout(1,1);
    set(gcf, 'Color',[1 1 1])
    ax1 = nexttile;
    imagesc(x_data,y_data,epidemic_size*100);
    if flip_yaxis_flag == true
        set(gca,'YDir','normal'); 
    end
    set(gca,'FontSize',tick_label_fontsize); 
    title(['Final size, ',variant]);
    yticks(yticks_vals); yticklabels(yticks_labels);
    c = colorbar; c.Ruler.TickLabelFormat='%g%%'; caxis([0,100]);
    xticks(xticks_vals); xticklabels(xticks_labels);
    ylabel(ylabel_text,'FontSize',label_fontsize);
    if ~isempty(xlabel_text)
        xlabel(xlabel_text,'FontSize',label_fontsize);
    end
    % if glyph_flag is active, add glyphs in designated positions
    if glyph_flag == true
        add_glyph_fn(glyph_pos,glyph_markertypes,glyph_colour,glyph_size,xticks_vals,yticks_vals)
    end
    ytickformat('%.2f')
    xtickformat('%.1f')
    colormap(gca,CT_purple)
    if save_figs_flag == true
        export_fig([save_filename_prefix,'_final_size'],'-pdf','-r1200','-transparent')
    end

    % Peak height
    f = figure;
    f.Position = [ 3   40   560*1.2   420*1.5];
    t = tiledlayout(1,1);
    set(gcf, 'Color',[1 1 1])
    ax1 = nexttile;
    imagesc(x_data,y_data,peak_height*100); 
    if flip_yaxis_flag == true
        set(gca,'YDir','normal'); 
    end
    set(gca,'FontSize',tick_label_fontsize); 
    title(['Peak infectious prevalence, ',variant]);
    xticks(xticks_vals); xticklabels(xticks_labels);  
    c = colorbar;c.Ruler.TickLabelFormat='%g%%'; caxis([0,10]);
    if ~isempty(xlabel_text)
        xlabel(xlabel_text,'FontSize',label_fontsize);
    end
    yticks(yticks_vals); yticklabels(yticks_labels);
    ylabel(ylabel_text,'FontSize',label_fontsize);
    % if glyph_flag is active, add glyphs in designated positions
    if glyph_flag == true
        add_glyph_fn(glyph_pos,glyph_markertypes,glyph_colour,glyph_size,xticks_vals,yticks_vals)
    end
    ytickformat('%.2f')
    xtickformat('%.1f')
    colormap(gca,CT_purple)
    box on
    if save_figs_flag == true
        export_fig([save_filename_prefix,'_peak_height'],'-pdf','-r1200','-transparent')
    end

    % Peak timing
    f = figure;
    f.Position = [ 3   40   560*1.2   420*1.5];
    t = tiledlayout(1,1);
    set(gcf, 'Color',[1 1 1])
    ax1 = nexttile;
    imagesc(x_data,y_data,datenum(peak_time)-datenum(date1)+1); 
    if flip_yaxis_flag == true
        set(gca,'YDir','normal'); 
    end
    set(gca,'FontSize',tick_label_fontsize); 
    title(['Time of peak prevalence, ',variant]);
    xticks(xticks_vals); xticklabels(xticks_labels); 
    cb = colorbar; formatOut = 'dd/mm/yy';  caxis([1,366]);
    cb.TickLabels = datestr(datetime(cb.Ticks-1+date1,'ConvertFrom','datenum'),formatOut);
    yticks(yticks_vals); yticklabels(yticks_labels);
    ylabel(ylabel_text,'FontSize',label_fontsize);
    if ~isempty(xlabel_text)
        xlabel(xlabel_text,'FontSize',label_fontsize);
    end
    % if glyph_flag is active, add glyphs in designated positions
    if glyph_flag == true
        add_glyph_fn(glyph_pos,glyph_markertypes,glyph_colour,glyph_size,xticks_vals,yticks_vals)
    end
    ytickformat('%.2f')
    xtickformat('%.1f')
    colormap(gca,CT_purple)
    box on
    if save_figs_flag == true
        export_fig([save_filename_prefix,'_peak_timing'],'-pdf','-r1200','-transparent')
    end
end

% Function to add glyph markers to heatmaps in Figure 1
function add_glyph_fn(glyph_pos,glyph_markertypes,glyph_colour,glyph_size,xticks_vals,yticks_vals)
    
    % Get spacing between each xtick and ytick data entry
    x_resolution = xticks_vals(2) - xticks_vals(1);
    
    % Set offset for text placement
    x_offset = [x_resolution/3 x_resolution/6 (x_resolution*(3/5) - 0.012)];
    y_offset = 0;

    % Default is the four corners of the heatmap
    % If glyph_pos is non-empty, plot in other positions
    if isempty(glyph_pos)
        % Put glpyhs in four corners
        text(xticks_vals(1),yticks_vals(1),glyph_markertypes{1},'Color',glyph_colour,'FontSize',glyph_size)
        text(xticks_vals(end),yticks_vals(1),glyph_markertypes{2},'Color',glyph_colour,'FontSize',glyph_size)
        text(xticks_vals(1),yticks_vals(end),glyph_markertypes{3},'Color',glyph_colour,'FontSize',glyph_size)
        text(xticks_vals(end),yticks_vals(end),glyph_markertypes{4},'Color',glyph_colour,'FontSize',glyph_size)
    else
        for glyph_itr = 1:size(glyph_pos,1)
            text(xticks_vals(glyph_pos(glyph_itr,1))-x_offset(glyph_itr),...
                 yticks_vals(glyph_pos(glyph_itr,2))-y_offset,...
                 glyph_markertypes{glyph_itr},...
                 'Color','r',...
                 'FontSize',glyph_size)
        end
    end
end

%% Function to plot Figure 3 heatmaps
function make_heatmap_plot_new_vacc(data,...
                                    x_vals,...
                                    y_vals,...
                                    label_fontsize,...
                                    xlabel_string,...
                                    ylabel_string,...
                                    title_string,...
                                    flip_yaxis_flag,...
                                    y_vals_is_date_type,...
                                    cbar_type,...
                                    cbar_range,...
                                    CT_map,...
                                    save_filename)

    % Set up figure
    f = figure;
    f.Position = [ 3   40   560*1.8   420*1.8];
    t = tiledlayout(1,1);
    set(gcf, 'Color',[1 1 1])

    % Generate heatmap
    data(isnan(data)) = -1; % Set NaNs as negative
    imagesc(x_vals,y_vals,data);

    % Flip y-axis if applicable
    if flip_yaxis_flag == true
        set(gca,'YDir','normal');
    end

    % Set title
    title(title_string);

    % If plotting peak time data, alter values to datetime format
    if y_vals_is_date_type == true
        formatOut_date_yaxis = 'ddmmmyyyy';  
        y_vals_labels = datestr(datetime(y_vals,'ConvertFrom','datenum'),formatOut_date_yaxis);
    else
        % No amendment to y_vals
        y_vals_labels = y_vals;
    end

    % Set axes labels
    yticks(y_vals); yticklabels(y_vals_labels);
    xticks(x_vals); xticklabels(x_vals);
    ylabel(ylabel_string,'FontSize',label_fontsize);
    xlabel(xlabel_string,'FontSize',label_fontsize);

    % X-axis: Set tick format and alter tick frequency
    xtickformat('%.2f')
    g = gca;
    if y_vals_is_date_type == true
        g.XTick = g.XTick(1:2:end);
    else
        % Epidemic probability plot
        g.XTick = g.XTick(2:2:end);
    end

    % Y-axis: Set tick format and alter tick frequency
    % Done if not date type labelling (epidemic probability plot)
    if y_vals_is_date_type == false
        ytickformat('%.1f')
    end

    % Add colourbar
    c = colorbar;
    if cbar_type == 1
        % Epidemic probability runs. Percentage that take off.
        c.Ruler.TickLabelFormat='%g%%';
        caxis(cbar_range);
    elseif cbar_type == 2
        % Time to attain specified prevalence
    %     caxis([0,365]);
    %     c.TickLabels(1) = {'NA'}; % Update description for scenario sets that had no outbreaks to NA
        formatOut = 'mmmyyyy';  
        %caxis([datenum(2021,5,17),datenum(2022,5,16)]);
        caxis(cbar_range);
        c.Ticks = [datenum(2021,6,1),datenum(2021,7,1),datenum(2021,8,1),...
                   datenum(2021,9,1),datenum(2021,10,1),datenum(2021,11,1),...
                   datenum(2021,12,1),datenum(2022,1,1),datenum(2022,2,1),...
                   datenum(2022,3,1),datenum(2022,4,1),datenum(2022,5,1)];
        c.TickLabels = datestr(datetime(c.Ticks,'ConvertFrom','datenum'),formatOut);
    elseif cbar_type==3
        % Percentages
        c.Ruler.TickLabelFormat='%g%%';
        caxis(cbar_range);
    elseif cbar_type==4
        % Epidemic probability
        c.Ruler.TickLabelFormat='%.2f';
        caxis(cbar_range);
    else
        error('Invalid cbar_type value.')
    end

    % Set colourmap
    CT_map(CT_map<0) = 0;
    % if cbar_type == 2
    %     CT_map(1,:) = [0.5 0.5 0.5];
    % end
    colormap(gca,CT_map)

    % Alter plot properties
    set(gca,'FontSize',label_fontsize);
    set(gca,'LineWidth',1);

    % Save figure
    if ~isempty(save_filename)
        export_fig(save_filename,'-pdf','-r1200')
    end
end

% Cloud plots/scatter plots of summary statistics from outbreak data
% Used for Figure S4
function generate_summary_measure_cloud_plots(epidemic_size_data,...
                                                peak_height_data,...
                                                peak_time_data,...
                                                colour_array,...
                                                markertypes,...
                                                leg_labels,...
                                                leg_pos,...
                                                fig_fontsize,...
                                                save_filename_resident_variants,...
                                                save_filename_VOC,...
                                                save_figs_flag)
    % Initialise the figures
    for fig_itr = 1:6
        position = [10, 10, 1.5*550, 1.5*450];
        set(0, 'DefaultFigurePosition', position);
        fig = figure(fig_itr);
        clf;
        set(fig,'Color', [1 1 1])
        hold on
    end

    % Set up the size of each marker
    sz = 50;

    % Get number of points being plotted
    n_data_points = size(epidemic_size_data,1);

    % Get colour maps to use on scatter plot
    n_colours = size(colour_array,1);
    CT_cell = cell(n_colours,1);
    for colour_itr = 1:n_colours

        % Get baseline colour RGB vector
        baseline_RGB = colour_array(colour_itr,:);

        % Create the colour transition colourmap
        % Go from white to the baseline_RGB
        CT_cell{colour_itr} = [linspace(baseline_RGB(1), 0.99, n_data_points)'... % Red colour component
                            linspace(baseline_RGB(2), 0.99, n_data_points)'... % Green colour component
                            linspace(baseline_RGB(3), 0.99, n_data_points)'];  % Blue colour component
    end

    % Populate the scatter plots
    for VOC_itr = 1:3

        %%% VOC scenarios %%%                
        figure(1)
        add_to_scatter_plot(epidemic_size_data(:,VOC_itr,2),...
                            peak_height_data(:,VOC_itr,2),...
                            sz,...
                            CT_cell{VOC_itr},...
                            markertypes{VOC_itr})

        figure(2)
        add_to_scatter_plot(peak_time_data(:,VOC_itr,2),...
                            epidemic_size_data(:,VOC_itr,2),...
                            sz,...
                            CT_cell{VOC_itr},...
                            markertypes{VOC_itr})

        figure(3)
        add_to_scatter_plot(peak_time_data(:,VOC_itr,2),...
                            peak_height_data(:,VOC_itr,2),...
                            sz,...
                            CT_cell{VOC_itr},...
                            markertypes{VOC_itr})

        %%% Resident variants %%%                
        figure(4)
        add_to_scatter_plot(epidemic_size_data(:,VOC_itr,1),...
                            peak_height_data(:,VOC_itr,1),...
                            sz,...
                            CT_cell{VOC_itr},...
                            markertypes{VOC_itr})

        figure(5)
        add_to_scatter_plot(peak_time_data(:,VOC_itr,1),...
                            epidemic_size_data(:,VOC_itr,1),...
                            sz,...
                            CT_cell{VOC_itr},...
                            markertypes{VOC_itr})

        figure(6)
        add_to_scatter_plot(peak_time_data(:,VOC_itr,1),...
                            peak_height_data(:,VOC_itr,1),...
                            sz,...
                            CT_cell{VOC_itr},...
                            markertypes{VOC_itr})
    end


    % For aligning axes limits, get limits for figures 1-3
    figure(1)
    xl_fig1 = xlim;
    yl_fig1 = ylim;

    figure(2)
    yl_fig2 = ylim;

    figure(3)
    yl_fig3 = ylim;

    % Set figure proerties
    for fig_itr = 1:6

        % Set axes labels & axes limits
        figure(fig_itr)
        if (fig_itr == 1) || (fig_itr == 4)
            xlabel('Outbreak size (propn of population)')
            ylabel('Peak infectious prevalence (propn of population)')
            xlim(xl_fig1)
            ylim(yl_fig1)
            ytickformat('%.3f')
        elseif (fig_itr == 2) || (fig_itr == 5)
            xlabel('Time of peak infectious prevalence')
            ylabel('Outbreak size (propn of population)')
            xlim([datetime(2021,5,1)  datetime(2022,6,1)])
            ylim(yl_fig2)
            xticks([datetime(2021,6,1) datetime(2021,9,1) datetime(2021,12,1) datetime(2022,3,1) datetime(2022,6,1)])
        elseif (fig_itr == 3) || (fig_itr == 6)
            xlabel('Time of peak infectious prevalence')
            ylabel('Peak infectious prevalence (propn of population)')
            xlim([datetime(2021,5,1)  datetime(2022,6,1)])
            xticks([datetime(2021,6,1) datetime(2021,9,1) datetime(2021,12,1) datetime(2022,3,1) datetime(2022,6,1)])
            ylim(yl_fig3)
            ytickformat('%.3f')
        end

        % Add title
        if fig_itr <= 3
            title('Variant of concern')
        else
            title('Resident variants')
        end

        % Add legend
        if (fig_itr == 3) || (fig_itr == 6)
            [~,icons] = legend(leg_labels,...
                    'Position',leg_pos,...
                    'FontSize',fig_fontsize-2);

            % Alter markersize in the legend
            for icon_itr = ((numel(icons)/2)+1):numel(icons)
                icons(icon_itr).Children.MarkerSize = 12;
            end
        end

        % Plot properties
        set(gca,'FontSize',fig_fontsize)
        set(gca,'LineWidth',1)
        box on

        % Save figure
        if save_figs_flag == true
            if fig_itr == 1
                export_fig([save_filename_VOC,'_outbreak_size_vs_peak_height'],'-pdf','-painters','-r1200')
            elseif fig_itr == 2
                export_fig([save_filename_VOC,'_outbreak_size_vs_peak_timing'],'-pdf','-painters','-r1200')
            elseif fig_itr == 3
                export_fig([save_filename_VOC,'_peak_timing_vs_peak_height'],'-pdf','-painters','-r1200')
            elseif fig_itr == 4
                export_fig([save_filename_resident_variants,'_outbreak_size_vs_peak_height'],'-pdf','-painters','-r1200')
            elseif fig_itr == 5
                export_fig([save_filename_resident_variants,'_outbreak_size_vs_peak_timing'],'-pdf','-painters','-r1200')
            elseif fig_itr == 6
                export_fig([save_filename_resident_variants,'_peak_timing_vs_peak_height'],'-pdf','-painters','-r1200')
            end
        end
    end
end
    % Sub-function used in scatter/cloud plot generation
function add_to_scatter_plot(x_data,y_data,sz,colour_vec,marker_type)
    if (marker_type == 's') || (marker_type == 'o')
        scatter(x_data,y_data,...
                sz,...
                colour_vec,...
                'filled',...
                'Marker',marker_type)
    else
        scatter(x_data,y_data,...
                sz,...
                colour_vec,...
                'Marker',marker_type)
    end    
end

% If you have the file_exchange function export_fig, then comment this out
function export_fig(savename,format,renderer,resolution)
exportgraphics(gcf,[savename,'.',format(2:end)])%,'resolution',resolution(2:end))
end

% Plot VOC scenario temporal proportion plots
function temporal_proportion_plots(outputs,...
                                    I_UK,...
                                    I_VOC,...
                                    line_type,...
                                    line_colours,...
                                    line_width,...
                                    glyph_markertypes,...
                                    xlim_end,...
                                    leg_flag,...
                                    leg_labels,...
                                    leg_pos,...
                                    save_figs_flag,...
                                    save_filename,...
                                    fig_fontsize,...
                                    vacc_data,...
                                    add_vaccination_data_flag)

    % Set x tick values to be used in upcoming plots
    xtick_vals = datetime(2021,6,1) + calmonths(0:2:11);

    % Plot proportions 
    position = [10, 10, 2.8*550, 1.5*450];
    set(0, 'DefaultFigurePosition', position);
    fig = figure(101);
    clf;
    set(fig,'Color', [1 1 1])
    hold on

    % Add to figure. Output depends on index of loop.
    y = I_VOC(:,:,1)*100./(I_VOC(:,:,1)+I_UK(:,:,1));
    x = outputs.dates;

    % Specify plotting against left y-axis
    yyaxis left
    ax_left = gca;

    % Add profiles to figure for each VOC
    for VOC_itr = 1:size(y,2)
        % set colour for line
        line_colour_vec = line_colours(VOC_itr,:);

        % Add line profile for VOC
        plot(x,...
                y(:,VOC_itr),...        
                line_type{VOC_itr},...
                'Color',line_colour_vec,...
                'LineWidth',line_width(VOC_itr),...
                'Marker',glyph_markertypes{VOC_itr},...
                'MarkerSize',12,...
                'MarkerIndices',1:20:length(y(:,VOC_itr)))
    end

    % y-axis properties
    ylim([0 105])

    % Set axis colour 
   ax_left.YAxis(1).Color = [0 0 0];

   % Panel specific
   ylabel('Percentage of cases that are VOC')

    % Set y-axis tick label style
    ytickformat('%.1f%%')

    % x-axis properties
    xlabel('Time');
    xlim([outputs.dates(1),xlim_end])
    xticks(xtick_vals)

    % Add vaccination data if required
    if add_vaccination_data_flag == true
        % Add vaccination onto the same panel
        add_vaccination_data_same_panel(outputs.dates,...
            vacc_data,...
            fig_fontsize,...
            true)
    end

    % Add legend
    if (leg_flag == true) %&& (fig_itr == 1)
        leg = legend(leg_labels,...
                        'Position',leg_pos);

        % Append vaccination bar data descriptor to legend
        bar_colour = rgba2rgb([1 1 1],[0 0 0.8 0.3]);
        bar(NaT,NaN,...
                        'FaceColor',bar_colour,...
                        'EdgeColor',bar_colour,...
                        'LineWidth',0.01,...
                        'DisplayName','Vaccination uptake (%)');

        % Set up indexing to move final label to be first, and shift all other
        % labels down by one
        n_labels = numel(leg.PlotChildren);
        label_reorder_idx = zeros(n_labels,1);
        label_reorder_idx(1) = n_labels;
        label_reorder_idx(2:end) = 1:(n_labels-1);

        % Reorder the legend labels
        reorder_leg = leg.PlotChildren(label_reorder_idx);
        legend(reorder_leg,'Position',leg_pos);
    end

    % Set plot properties
    set(ax_left,'Fontsize',fig_fontsize)
    set(ax_left,'LineWidth',1)
    box on

    % Save figures
    if save_figs_flag == true
        export_fig(save_filename,'-pdf','-painters','-r1200')
    end
end


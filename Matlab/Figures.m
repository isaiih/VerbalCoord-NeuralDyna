%% Brain representation - Global effect (VP) & correlation with PLV ; Phase-Amplitude Coupling %%


%% Initialisation 

addpath(['add path' '\fieldtrip-20201103']);
ft_defaults;

%% Electrodes localization %% (Figure 1)
    
path_load = 'add path here';
load([path_load 'Alls_hfa.mat'] );
load([path_load 'Allpop_hfa.mat'] );

% Split in left and right COIs
cfg = [];
cfg.channel =  '*''_*'; %'*''*-'; 
Allpop_L = ft_selectdata(cfg,Allpop);
cfg.channel = {'all'; '-*''_*'};
Allpop_R = ft_selectdata(cfg,Allpop);

for s = 1:length(idx_Pat)
    eval(['tmp = [AllS.' idx_Pat{s,1} '.elec];']);

    coi_L = []; coi_R = [];

    for i=1:length(tmp.label)
        if contains(tmp.label{i},'''')
            coi_L = [coi_L i];
        else
            coi_R = [coi_R i];
        end
    end
    tmp_L.coordsys = tmp.coordsys; tmp_L.unit = tmp.unit; tmp_L.label = tmp.label(coi_L);  tmp_L.elecpos = tmp.elecpos(coi_L,:);
    tmp_R.coordsys = tmp.coordsys; tmp_R.unit = tmp.unit; tmp_R.label = tmp.label(coi_R);  tmp_R.elecpos = tmp.elecpos(coi_R,:);

    eval(['AllS_L.' idx_Pat{s,1} '= tmp_L;']);
    eval(['AllS_R.' idx_Pat{s,1} '= tmp_R;']);

    clear tmp* coi*;
end

% get ft pial template
[ftver, path_ft] = ft_version;
load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;

% Set the colormap
clrmp = lines;
clrmp(8,1) = 0.75;  clrmp(9,1) = 0.95;  clrmp(10,1) = 1;
clrmp(8,2) = 0.15;  clrmp(9,2) = 0.95;  clrmp(10,2) = 0.1;
clrmp(8,3) = 0.92;  clrmp(9,3) = 0.95;  clrmp(10,3) = 0.1;

clrmp(11,1) = 0.1;  clrmp(12,1) = 0.1;  clrmp(13,1) = 0.1;
clrmp(11,2) = 1;    clrmp(12,2) = 0.15; clrmp(13,2) = 0.1; 
clrmp(11,3) = 0.1;  clrmp(12,3) = 1;    clrmp(13,3) = 0.1;

clrmp(14,1) = 1;    clrmp(15,1) = 1;    clrmp(16,1) = 0.86;
clrmp(14,2) = 1;    clrmp(15,2) = 0;    clrmp(16,2) = 0.48;
clrmp(14,3) = 0;    clrmp(15,3) = 1;    clrmp(16,3) = 0.15;

% Figure - side view
figure('Position', [50 50 1200 600]); hold on;
ax_1 = subplot(1,2,1);
ft_plot_mesh(fspial_lh,'facealpha', .10);
view([-120 0]); lighting gouraud; camlight HEADLIGHT;
L_Pat=[];
for s = 1:length(idx_Pat)
    % pick Left cortical
    eval(['tmp = [AllS_L.' idx_Pat{s,1} '];']);
    if ~isempty(tmp.elecpos)
        ft_plot_sens(tmp, 'facecolor', clrmp(s,:), 'elecshape', 'sphere', 'elecsize', 3,'edgecolor', 'k', 'edgealpha', .1);
        L_Pat = [L_Pat, s];
    end    
end
legend('elec hemi.left', idx_Pat{L_Pat,1}, 'Interpreter', 'none');
ax_2 = subplot(1,2,2);
ft_plot_mesh(fspial_rh,'facealpha', .10);
view([120 0]); lighting gouraud; camlight HEADLIGHT;
R_Pat=[];
for s = 1:length(idx_Pat)
    % pick Right cortical
    eval(['tmp = [AllS_R.' idx_Pat{s,1} '];']);
    if ~isempty(tmp.elecpos)
        ft_plot_sens(tmp, 'facecolor', clrmp(s,:), 'elecshape', 'sphere', 'elecsize', 3,'edgecolor', 'k', 'edgealpha', .1);
        R_Pat = [R_Pat, s]; 
    end
end
legend('elec hemi.right',idx_Pat{R_Pat,1}, 'Interpreter', 'none');
% path_sav = ;
% savefig(gcf,[path_sav '\elec_implantations_SIDE_VIEW'], 'compact');


%% Regions of interest - Visualisation %% (Figure 4 panel A)

% STG BA41/42 - STG BA22 - IPL BA40 - IFG BA44 %

roi    = {'STG, Left Superior Temporal Gyrus A41/42'; 'STG, Left Superior Temporal Gyrus A22'; ...  
          'IPL, Left Inferior Parietal Lobule A40'; 'IFG, Left Inferior Frontal Gyrus A44';};
clrmp = lines;
% get ft pial template
[ftver, path_ft] = ft_version;
load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;
% get atlas
% Brainnetome
atlas = ft_read_atlas([path_ft '/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii']);
atlas.coordsys = 'acpc';

for i=1:length(roi)
    % define left and right ROIs
    isLeft(i) = contains(roi{i}, 'left', 'IgnoreCase',true);
    isRight(i) = contains(roi{i}, 'right', 'IgnoreCase',true);
    if isLeft(i)==0 && isRight(i)==0
        error(['CHECK YOUR cfg' newlline 'roi name must specify the lateralization (left or right)']);
    end

    % find the parcel(s) corresponding to ROI input name
    roi_atlas = contains(atlas.tissuelabel, roi{i});
    roi_atlas = atlas.tissuelabel(roi_atlas);
    if length(roi_atlas)>1
        warning(['input ROI was: "' roi{i} '"' newline 'which correspond to ' num2str(length(roi_atlas)) ' parcel names:' newline ...
            roi_atlas{:}]);
    end

    % create a volumetric mask of the ROIs
    cfg            = [];
    cfg.inputcoord = 'acpc';
    cfg.atlas      = atlas;
    cfg.roi        = roi_atlas;
    mask           = ft_volumelookup(cfg, atlas);
    ROI         = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    ROI.brain   = mask;

    % create triangulated and smoothed surface mesh
    cfg             = [];
    cfg.method      = 'iso2mesh';
    cfg.radbound    = 2;
    cfg.maxsurf     = 0;
    cfg.tissue      = 'brain';
    cfg.numvertices = 100000;
    cfg.smooth      = 3;
    cfg.spmversion  = 'spm12';
    mesh{i}         = ft_prepare_mesh(cfg, ROI);
end

for i=1:length(roi)
    if isLeft(i)
        ft_plot_mesh(mesh{i},'facealpha', .209, 'facecolor', clrmp(i,:), 'edgecolor', 'none');
    end
end


%% Figures - Global effect (VP) %% (Figure 3 panel A ; Figure S2 panel A)

freq_names = {'delta'; 'theta'; 'alpha'; 'beta'; 'low gamma'; 'hfa'};

for choice_freq=1:length(freq_names)

    freq = freq_names{choice_freq};
    
    path_load = 'add path here';
    load([path_load 'Alls_' freq '.mat'] );
    load([path_load 'Allpop_' freq '.mat'] );
    load([path_load 'Allpop_L_' freq '.mat'] );
    load([path_load 'Allpop_R_' freq '.mat'] );

    % get ft pial template
    [ftver, path_ft] = ft_version;
    load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
    load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;

    % parameters for the colormap
    clrmp = jet;
    
    % to select only the significant electrodes (permutation test)
    % both hemispheres
    for i=1:height(Allpop.seeg)
        if Allpop.t_permut(i) == 0
            Allpop.seeg(i) = 0;
        end
    end
    % left hemisphere
    for i=1:height(Allpop_L.seeg)
        if Allpop_L.t_permut(i) == 0
            Allpop_L.seeg(i) = 0;
        end
    end    
    % right hemisphere
    for i=1:height(Allpop_R.seeg)
        if Allpop_R.t_permut(i) == 0
            Allpop_R.seeg(i) = 0;
        end
    end
    
    % to remove no significant channels
    % both hemispheres
    my_no_roi = {};
    poc = 1;
    for i=1:height(Allpop.seeg)
        if Allpop.seeg(i) == 0
            my_no_roi(poc,1) = Allpop.label(i);
            poc = poc + 1;
        end
    end
    for z=1:height(my_no_roi)
        host_len = height(Allpop.label) - 1;
        for i=1:host_len
            if Allpop.seeg(i) == 0 
                Allpop.label(i) = [];
                Allpop.seeg(i) = [];
                Allpop.t_permut(i) = [];
                break
            end
        end
    end
    if Allpop.seeg(end) == 0 
        Allpop.label(end) = [];
        Allpop.seeg(end) = [];
        Allpop.t_permut(end) = [];
    end
    %left
    my_no_roi_L = {};
    poc = 1;
    for i=1:height(Allpop_L.seeg)
        if Allpop_L.seeg(i) == 0
            my_no_roi_L(poc,1) = Allpop_L.label(i);
            poc = poc + 1;
        end
    end
    for z=1:height(my_no_roi_L)
        host_len = height(Allpop_L.label) - 1;
        for i=1:host_len
            if Allpop_L.seeg(i) == 0 
                Allpop_L.label(i) = [];
                Allpop_L.seeg(i) = [];
                Allpop_L.t_permut(i) = [];
                break
            end
        end
    end
    if Allpop_L.seeg(end) == 0 
        Allpop_L.label(end) = [];
        Allpop_L.seeg(end) = [];
        Allpop_L.t_permut(end) = [];
    end     
    %right
    my_no_roi_R = {};
    poc = 1;
    for i=1:height(Allpop_R.seeg)
        if Allpop_R.seeg(i) == 0
            my_no_roi_R(poc,1) = Allpop_R.label(i);
            poc = poc + 1;
        end
    end
    for z=1:height(my_no_roi_R)
        host_len = height(Allpop_R.label) - 1;
        for i=1:host_len
            if Allpop_R.seeg(i) == 0 
                Allpop_R.label(i) = [];
                Allpop_R.seeg(i) = [];
                Allpop_R.t_permut(i) = [];
                break
            end
        end
    end
    if Allpop_R.seeg(end) == 0 
        Allpop_R.label(end) = [];
        Allpop_R.seeg(end) = [];
        Allpop_R.t_permut(end) = [];
    end

    % Plot 
    % Both hemispheres
    cfg              = [];
    cfg.funcolormap  = clrmp;
    cfg.funparameter = 'seeg';
    cfg.method       = 'cloud';
    cfg.cloudtype    = 'point';
    cfg.scalerad     = 'no';
    cfg.radius       = 40;
    cfg.facealpha    = .1;
    cfg.frequency    = Allpop.freq(1); 
    cfg.funcolorlim  = [-1 1];
    ft_sourceplot(cfg, Allpop);
    ft_plot_mesh(fspial_lh,'facealpha', .15);
    ft_plot_mesh(fspial_rh,'facealpha', .15);
    view([90 0]); lighting gouraud; camlight('headlight', 'infinite');
    view([180 30]); lighting gouraud; camlight('headlight', 'infinite');
    view([-180 -90]); lighting gouraud; camlight('headlight', 'infinite');
    % path_sav = [];
    % savefig(gcf,[path_sav '\VP_global_effect_' freq '_BOTH'], 'compact');

    %left hemisphere
    ft_sourceplot(cfg, Allpop_L);
    ft_plot_mesh(fspial_lh,'facealpha', .10);
    view([-180 -90]); lighting gouraud; camlight('headlight', 'infinite');
    % path_sav = [];
    % savefig(gcf,[path_sav '\VP_global_effect_' freq '_LEFT'], 'compact');
    
    %right hemisphere
    ft_sourceplot(cfg, Allpop_R);
    ft_plot_mesh(fspial_rh,'facealpha', .10);
    view([180 -90]); lighting gouraud; camlight('headlight', 'infinite');
    % path_sav = [];
    % savefig(gcf,[path_sav '\VP_global_effect_' freq '_RIGHT'], 'compact');    
end


%% Figures - Correlation (VP global effect ~ synchr.speech PLV) %% (Figure 3 panel B ; Figure S2 panel B)

freq_names = {'delta'; 'theta'; 'alpha'; 'beta'; 'low gamma'; 'hfa'};

for choice_freq=1:length(freq_names)

    freq = freq_names{choice_freq};
    
    path_load = 'add path here';
    load([path_load 'Alls_' freq '.mat'] );
    load([path_load 'Allpop_' freq '.mat'] );
    load([path_load 'Allpop_L_' freq '.mat'] );
    load([path_load 'Allpop_R_' freq '.mat'] );

    % get ft pial template
    [ftver, path_ft] = ft_version;
    load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
    load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;


    % parameters for the colormap
    clrmp = jet;
    extend_blue_2 = transpose(linspace(0,1,97));
    clrmp(32:128,2) = extend_blue_2;
    clrmp(97:128,1) = 0;
    clrmp(97:128,3) = 1;
    extend_orange_2 = flip(transpose(linspace(0,1,97)));
    clrmp(128:224,2) = extend_orange_2;
    clrmp(128:159,1) = 1;
    clrmp(129:159,3) = 0;
      
    % to select the significant channels
    % both hemispheres
    for i=1:length(Allpop.t_permut_spearman) 
        if Allpop.t_permut_spearman(i) == 0  
            Allpop.rho_spearman(i) = 0;      
        end 
    end
    % left hemisphere
    for i=1:length(Allpop_L.t_permut_spearman) 
        if Allpop_L.t_permut_spearman(i) == 0  
            Allpop_L.rho_spearman(i) = 0;      
        end 
    end
    % right hemisphere
    for i=1:length(Allpop_R.t_permut_spearman) 
        if Allpop_R.t_permut_spearman(i) == 0  
            Allpop_R.rho_spearman(i) = 0;      
        end 
    end
    
    % to remove no significant channels
    % both hemispheres
    my_no_roi = {};
    poc = 1;
    for i=1:height(Allpop.rho_spearman)
        if Allpop.rho_spearman(i) == 0
            my_no_roi(poc,1) = Allpop.label(i);
            poc = poc + 1;
        end
    end
    for z=1:height(my_no_roi)
        host_len = height(Allpop.label) - 1;
        for i=1:host_len
            if Allpop.rho_spearman(i) == 0 
                Allpop.label(i) = [];
                Allpop.rho_spearman(i) = [];
                Allpop.t_permut_spearman(i) = [];
                break
            end
        end
    end
    if Allpop.rho_spearman(end) == 0 
        Allpop.label(end) = [];
        Allpop.rho_spearman(end) = [];
        Allpop.t_permut_spearman(end) = [];
    end
    % left hemisphere
    my_no_roi_L = {};
    poc = 1;
    for i=1:height(Allpop_L.rho_spearman)
        if Allpop_L.rho_spearman(i) == 0
            my_no_roi_L(poc,1) = Allpop_L.label(i);
            poc = poc + 1;
        end
    end
    for z=1:height(my_no_roi_L)
        host_len = height(Allpop_L.label) - 1;
        for i=1:host_len
            if Allpop_L.rho_spearman(i) == 0 
                Allpop_L.label(i) = [];
                Allpop_L.rho_spearman(i) = [];
                Allpop_L.t_permut_spearman(i) = [];
                break
            end
        end
    end
    if Allpop_L.rho_spearman(end) == 0 
        Allpop_L.label(end) = [];
        Allpop_L.rho_spearman(end) = [];
        Allpop_L.t_permut_spearman(end) = [];
    end    
    % right hemisphere
    my_no_roi_R = {};
    poc = 1;
    for i=1:height(Allpop_R.rho_spearman)
        if Allpop_R.rho_spearman(i) == 0
            my_no_roi_R(poc,1) = Allpop_R.label(i);
            poc = poc + 1;
        end
    end
    for z=1:height(my_no_roi_R)
        host_len = height(Allpop_R.label) - 1;
        for i=1:host_len
            if Allpop_R.rho_spearman(i) == 0 
                Allpop_R.label(i) = [];
                Allpop_R.rho_spearman(i) = [];
                Allpop_R.t_permut_spearman(i) = [];
                break
            end
        end
    end
    if Allpop_R.rho_spearman(end) == 0 
        Allpop_R.label(end) = [];
        Allpop_R.rho_spearman(end) = [];
        Allpop_R.t_permut_spearman(end) = [];
    end

    
    % Plot
    
    % both hemispheres
    cfg              = [];
    cfg.funcolormap  = clrmp;
    cfg.funparameter = 'rho_spearman';
    cfg.method       = 'cloud';
    cfg.cloudtype    = 'point';
    cfg.scalerad     = 'no';
    cfg.radius       = 40;
    cfg.facealpha    = .1;
    cfg.frequency    = Allpop.freq(1); 
    cfg.funcolorlim  = [-1 1];

    ft_sourceplot(cfg, Allpop);
    ft_plot_mesh(fspial_lh,'facealpha', .12);
    ft_plot_mesh(fspial_rh,'facealpha', .12);
    view([90 0]); lighting gouraud; camlight('headlight', 'infinite');
    view([-90 0]); lighting gouraud; camlight('headlight', 'infinite');
    view([-180 -90]); lighting gouraud; camlight('headlight', 'infinite');    
    view([180 30]); lighting gouraud; camlight('headlight', 'infinite');
    %Create Colorbar
    cbh = colorbar ; 
    cbh.Ticks = [-1,-0.5, -0.2, 0, 0.2, 0.5, 1]; 
    cbh.TickLabels = [-1, -0.5, -0.2, 0, 0.2, 0.5, 1];  
    cbh.FontSize = 12; 
    cbh.Label.String = "r values (power ~ PLV)";
    cbh.Label.FontSize = 15;
    cbh.Label.Position = [-2 0 0];
    % path_sav = [];
    % savefig(gcf,[path_sav 'PLV_correlation_' freq '_BOTH'], 'compact');

    % left hemisphere
    ft_sourceplot(cfg, Allpop_L);
    ft_plot_mesh(fspial_lh,'facealpha', .12);
    view([-90 0]); lighting gouraud; camlight('headlight', 'infinite');
    %Create Colorbar
    cbh = colorbar ; 
    cbh.Ticks = [-1,-0.5, -0.2, 0, 0.2, 0.5, 1]; 
    cbh.TickLabels = [-1, -0.5, -0.2, 0, 0.2, 0.5, 1];  
    cbh.FontSize = 12;
    cbh.Label.String = "r values (power ~ PLV)";
    cbh.Label.FontSize = 15;
    cbh.Label.Position = [-2 0 0];
    % path_sav = [];
    % savefig(gcf,[path_sav 'PLV_correlation_' freq '_LEFT'], 'compact');
    
    % right hemisphere
    ft_sourceplot(cfg, Allpop_R);
    ft_plot_mesh(fspial_rh,'facealpha', .12);
    view([90 0]); lighting gouraud; camlight('headlight', 'infinite');
    %Create Colorbar
    cbh = colorbar ; 
    cbh.Ticks = [-1,-0.5, -0.2, 0, 0.2, 0.5, 1]; 
    cbh.TickLabels = [-1, -0.5, -0.2, 0, 0.2, 0.5, 1];  
    cbh.FontSize = 12;
    cbh.Label.String = "r values (power ~ PLV)";
    cbh.Label.FontSize = 15;
    cbh.Label.Position = [-2 0 0];
    % path_sav = [];
    % savefig(gcf,[path_sav 'PLV_correlation_' freq '_RIGHT'], 'compact');    
end

%% Figures - PAC analyses %% (Figure 5 panel A)

% PAC - Virtual Partner speech (phase)

freq = 'hfa'; % used for the amplitude

path_load = 'add path here';
load([path_load 'Allpop_L_' freq '.mat'] );

% get ft pial template
[ftver, path_ft] = ft_version;
load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;

% parameters for the colormap
clrmp = flip(hot);
clrmp(1,:) = 1;

% Plot
% Left hemisphere
cfg              = [];
cfg.funcolormap  = clrmp;
cfg.funparameter = 'VP_pac_val';
cfg.method       = 'cloud';
cfg.cloudtype    = 'point';
cfg.scalerad     = 'no';
cfg.radius       = 30;
cfg.facealpha    = .5;
cfg.frequency    = Allpop_L.freq(1); 
cfg.funcolorlim  = [0.2  0.7];

ft_sourceplot(cfg, Allpop_L);
ft_plot_mesh(fspial_lh,'facealpha', .12);
view([-90 0]); lighting gouraud; camlight('headlight', 'infinite');
cbh = colorbar ; 
cbh.Ticks = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]; 
cbh.TickLabels = [20, 30, 40, 50, 60, 70];  
cbh.FontSize = 12;
cbh.Label.String = "significant PAC (VP speech)";
cbh.Label.FontSize = 15;
cbh.Label.Position = [-2 0 0];

% To add the region of interest - Inferior frontal Gyrus BA44 %
roi    = {'IFG, Left Inferior Frontal Gyrus A44';}; 
clrmp = lines; 
% get ft pial template
[ftver, path_ft] = ft_version;
load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;
% get atlas
% Brainnetome
atlas = ft_read_atlas([path_ft '/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii']);
atlas.coordsys = 'acpc';

for i=1:length(roi)

    % define left and right ROIs
    isLeft(i) = contains(roi{i}, 'left', 'IgnoreCase',true);
    isRight(i) = contains(roi{i}, 'right', 'IgnoreCase',true);
    if isLeft(i)==0 && isRight(i)==0
        error(['CHECK YOUR cfg' newlline 'roi name must specify the lateralization (left or right)']);
    end

    % find the parcel(s) corresponding to ROI input name
    roi_atlas = contains(atlas.tissuelabel, roi{i});
    roi_atlas = atlas.tissuelabel(roi_atlas);
    if length(roi_atlas)>1
        warning(['input ROI was: "' roi{i} '"' newline 'which correspond to ' num2str(length(roi_atlas)) ' parcel names:' newline ...
            roi_atlas{:}]);
    end

    % create a volumetric mask of the ROIs
    cfg            = [];
    cfg.inputcoord = 'acpc';
    cfg.atlas      = atlas;
    cfg.roi        = roi_atlas;
    mask           = ft_volumelookup(cfg, atlas);
    ROI         = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    ROI.brain   = mask;

    % create triangulated and smoothed surface mesh
    cfg             = [];
    cfg.method      = 'iso2mesh';
    cfg.radbound    = 2;
    cfg.maxsurf     = 0;
    cfg.tissue      = 'brain';
    cfg.numvertices = 100000;
    cfg.smooth      = 3;
    cfg.spmversion  = 'spm12';
    mesh{i}         = ft_prepare_mesh(cfg, ROI);

end

ft_plot_mesh(fspial_lh, 'facealpha', .12, 'facecolor', 'skin', 'edgecolor', 'none');
for i=1:length(roi)
    if isLeft(i)
        ft_plot_mesh(mesh{i},'facealpha', .1, 'facecolor', clrmp(i,:), 'edgecolor', 'none');
    end
end
% path_sav = [];
% savefig(gcf,[path_sav 'VP_' freq '_LEFT'], 'compact');


% PAC - Phase difference (dynamic coordination/synchronisation)

freq = 'hfa'; % used for the amplitude

path_load = 'add path here';
load([path_load 'Allpop_L_' freq '.mat'] );

% get ft pial template
[ftver, path_ft] = ft_version;
load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;

% parameters for the colormap
clrmp = flip(hot);
clrmp(1,:) = 1;

% Plot
% Left hemisphere
cfg              = [];
cfg.funcolormap  = clrmp;
cfg.funparameter = 'PHDIFF_pac_val';
cfg.method       = 'cloud';
cfg.cloudtype    = 'point';
cfg.scalerad     = 'no';
cfg.radius       = 30;
cfg.facealpha    = .5;
cfg.frequency    = Allpop_L.freq(1); 
cfg.funcolorlim  = [0.2  0.7];

ft_sourceplot(cfg, Allpop_L);
ft_plot_mesh(fspial_lh,'facealpha', .12);
view([-90 0]); lighting gouraud; camlight('headlight', 'infinite');
cbh = colorbar ; 
cbh.Ticks = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]; 
cbh.TickLabels = [20, 30, 40, 50, 60, 70]; 
cbh.FontSize = 12;
cbh.Label.String = "significant PAC (Phase diff.)";
cbh.Label.FontSize = 15;
cbh.Label.Position = [-2 0 0];

% To add the region of interest - Inferior frontal Gyrus BA44 %
roi    = {'IFG, Left Inferior Frontal Gyrus A44';}; 
clrmp = lines; 
% get ft pial template
[ftver, path_ft] = ft_version;
load([path_ft filesep 'template\anatomy\surface_pial_left.mat']); fspial_lh = mesh;
load([path_ft filesep 'template\anatomy\surface_pial_right.mat']); fspial_rh = mesh; clear mesh;
% get atlas
% Brainnetome
atlas = ft_read_atlas([path_ft '/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii']);
atlas.coordsys = 'acpc';

for i=1:length(roi)

    % define left and right ROIs
    isLeft(i) = contains(roi{i}, 'left', 'IgnoreCase',true);
    isRight(i) = contains(roi{i}, 'right', 'IgnoreCase',true);
    if isLeft(i)==0 && isRight(i)==0
        error(['CHECK YOUR cfg' newlline 'roi name must specify the lateralization (left or right)']);
    end

    % find the parcel(s) corresponding to ROI input name
    roi_atlas = contains(atlas.tissuelabel, roi{i});
    roi_atlas = atlas.tissuelabel(roi_atlas);
    if length(roi_atlas)>1
        warning(['input ROI was: "' roi{i} '"' newline 'which correspond to ' num2str(length(roi_atlas)) ' parcel names:' newline ...
            roi_atlas{:}]);
    end

    % create a volumetric mask of the ROIs
    cfg            = [];
    cfg.inputcoord = 'acpc';
    cfg.atlas      = atlas;
    cfg.roi        = roi_atlas;
    mask           = ft_volumelookup(cfg, atlas);
    ROI         = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    ROI.brain   = mask;

    % create triangulated and smoothed surface mesh
    cfg             = [];
    cfg.method      = 'iso2mesh';
    cfg.radbound    = 2;
    cfg.maxsurf     = 0;
    cfg.tissue      = 'brain';
    cfg.numvertices = 100000;
    cfg.smooth      = 3;
    cfg.spmversion  = 'spm12';
    mesh{i}         = ft_prepare_mesh(cfg, ROI);

end

for i=1:length(roi)
    if isLeft(i)
        ft_plot_mesh(mesh{i},'facealpha', .1099, 'facecolor', clrmp(i,:), 'edgecolor', 'none');
    end
end
% path_sav = [];
% savefig(gcf,[path_sav 'PHDIFF_' freq '_LEFT'], 'compact');



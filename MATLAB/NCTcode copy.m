%% SETUP

% Specify the folder containing the .mat files
folderPath = '/Users/daraneumann/Downloads/Research/NCT';
addpath '/Users/daraneumann/Downloads/Research/NCT'

% Get a list of all .mat files in the folder
matFiles = dir(fullfile(folderPath, '*.mat'));
newTS = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilTimeSeriesZScored.mat').labpcexport;

%% COMPUTE BEST NUMBER OF CLUSTERS

% Set file location and save directory
basedir = '/Users/daraneumann/Downloads/Research/NCT/';
cd(basedir); addpath(genpath('code'))
savedir = fullfile(basedir,'Results','Clusters');mkdir(savedir);

% Load the Z-scored time series data
% TS = load('/Users/daraneumann/Downloads/Research/NCT/Data/TimeSeriesZScored.mat').labpcexport;
% size(TS)

% Set parameters
distanceMethod = 'correlation'; % Distance metric
nreps = 50; % Number of times to repeat clustering; lowest error solution will be chosen
maxI = 1000; % Tolerance for k-means convergence attempts
split = 10; % Split amount


% Run the loop (runs for approximately 5 hours)
for numClusters = 2:10
    
    disp(['Starting K = ', num2str(numClusters), ' Split = ',num2str(split)]);
    
    [partition,~,sumd] = kmeans(newTS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    save(fullfile(savedir,['KMeans.TimeSeriesZScored',num2str(split),'.K',num2str(numClusters),'.mat']),'partition','sumd')
    
    disp(['Finished K = ',num2str(numClusters),' Split = ',num2str(split)]);
    
    clear partition
    clear sumd
    
    
end

disp('Done computing clusters.');

%% VARIANCE EXPLAINED BY CLUSTERS

% Locate the MATLAB function in the EnergyLandscape folder
addpath /Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/assesscluster/

% Initialize values for the loop
k_range = 2:10;
VarianceExplained = zeros(length(k_range),1);

% Run the loop
for numClusters = k_range
    disp(['K = ',num2str(numClusters)])
    load(fullfile(basedir,['/Results/Clusters/KMeans.TimeSeriesZScored',num2str(split),'.K',num2str(numClusters),'.mat']));
    
    kClusterCentroids = GET_CENTROIDS(newTS,partition,numClusters);
    VarianceExplained(numClusters - 1) = VAREXPLAINED(newTS,partition,kClusterCentroids,numClusters);
end

% Save one of the variance explained files (10 clusters in this case)
save(fullfile(savedir,['AprilVarianceExplained.kmeans.TimeSeriesZScored',num2str(split),'.mat']),'VarianceExplained','k_range');


% Create figures
f=figure;
plot(k_range,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
set(f, 'PaperUnits', 'inches');
x_width=6 ;y_width=6;
set(f, 'PaperPosition', [0 0 x_width y_width]);
saveas(f,'/Users/daraneumann/Downloads/Research/NCT/Results/Figures/AprilVariance.Explained.Graph.kmeans.TS.ZScored.png');

f = figure;
plot(k_range(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
set(f, 'PaperUnits', 'inches');
x_width=6 ;y_width=6;
set(f, 'PaperPosition', [0 0 x_width y_width]); %
saveas(f,'/Users/daraneumann/Downloads/Research/NCT/Results/Figures/AprilGain.Variance.Explained.Graph.kmeans.TS.ZScored.png');


%% PARTITIONS AND CENTROIDS

% Locate MATLAB function
addpath /Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/assesscluster/

% Set parameters
numClusters = 6; % Best number of clusters
distanceMethod = 'correlation'; % Distance metric
nreps = 50; % Number of times to repeat clustering (same as previous)

% Compute partition (runs for approximately 15 minutes)
[Partition_ZS_TS,~,sumd] = kmeans(newTS,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',1000);
size(Partition_ZS_TS)
Partition_k6_ZS_TS = Partition_ZS_TS;

% Save and load partition
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilPartition.k6.ZScored.TS.mat','Partition_ZS_TS');
% Partition_k6_ZS_TS = load('/Users/daraneumann/Downloads/Research/NCT/Results/Partition.k6.ZScored.TS.mat').Partition_ZS_TS;


% Compute centroids
Centroids_k6_ZS_TS = GET_CENTROIDS(newTS,Partition_ZS_TS,numClusters);


% Save and load centroid
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilCentroids.k6.ZScored.TS.mat','Centroids_k6_ZS_TS');
% Centroids_k6_ZS_TS = load('/Users/daraneumann/Downloads/Research/NCT/Results/Centroids.k6.ZScored.TS.mat').Centroids_k6_ZS_TS;


%% TIME SERIES LENGTH

% Load partition and centroids
% Partition_k6_ZS_TS = load('/Users/daraneumann/Downloads/Research/NCT/Results/AprilPartition.k6.ZScored.TS.mat').Partition_ZS_TS;
% Centroids_k6_ZS_TS = load('/Users/daraneumann/Downloads/Research/NCT/Results/Centroids.k6.ZScored.TS.mat').Centroids_k6_ZS_TS;

% Import time series length data
length_timeseries = load('/Users/daraneumann/Downloads/Research/NCT/Data/LengthTS.mat');
length_timeseries = length_timeseries.labpcexport;

% Identify best number of clusters
best_number_of_clusters = max(Partition_k6_ZS_TS);
disp(best_number_of_clusters)

% Assign subject IDs
n = 554;
subjID = [];
for i = 1:n
    leng = length_timeseries(i);
    subjID = [subjID; ones(leng, 1)*i]; %#ok<AGROW>
end


%% TRANSITION PROBABILITIES

% Locate the MATLAB function
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/dynamicsfxns'

% Initialize values
best_number_of_clusters = 6;
n = 554;
[transProbs,transitionProbabilityMats,numTransitions] = GET_TRANS_PROBS(Partition_k6_ZS_TS, subjID, best_number_of_clusters);
tran_prob = cell(n, 1);

for i = 1:n
    vec = i;
    tmp = transitionProbabilityMats(vec,:,:);
    for j = 1
        for p = 1:best_number_of_clusters
            tran_prob{i,j} = squeeze(tmp(j,:,:));
        end
    end
end

% Save transition probability
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilTran.Prob.k6.ZS.TS.mat','tran_prob')


%% CENTROIDS FOR ALL PATIENTS


addpath /Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/assesscluster/


% Initialize values
numClusters = 6;
n = 554;

% Run a loop to compile all centroids
Centroids_PerPt_ZS_TS = zeros(n,200,numClusters);

for subjectno = 1:n
    row = find(subjID==subjectno);
    Centroids_PerPt_ZS_TS(subjectno,:,:) = GET_CENTROIDS(newTS(row,:),Partition_k6_ZS_TS(row),numClusters);
end

size(Centroids_PerPt_ZS_TS);

% Validate size of the matrix
size(Centroids_PerPt_ZS_TS)

% Validate SC map of sample subjects
% imagesc(reshape(Centroids_PerPt_ZS_TS(1,:,:),200,numClusters));

% Save and load values of all centroids
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilCentroids.PerPt.ZS.TS.mat','Centroids_PerPt_ZS_TS')
% Centroids_PerPt_ZS_TS=load('/Users/daraneumann/Downloads/Research/NCT/Results/Centroids.PerPt.ZS.TS.mat').Centroids_PerPt_ZS_TS;

% Create an SC map of the averaged centroids across all subjects
Averaged_CentroidsZS_TS = reshape(mean(Centroids_PerPt_ZS_TS(1:n,:,:)),200,numClusters);imagesc(Averaged_CentroidsZS_TS);

%% FIND BEST T

% Initialize values
nperms = 1000;
SC = load('/Users/daraneumann/Downloads/Research/NCT/Data/SC.Matrix.diag0.mat').labpcexport;
numClusters = 6;

% Locate the MATLAB function
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/control'

% Set parameters
c = 1; % Time scale parameter based on values from paper
Anorm = NORMALIZE(SC,c); % Normalize A by maximum eigenvalue - eye(N) to make marginally stable
T_range = 0.001:0.5:10; nT = length(T_range);
Xf_ind = repmat(1:numClusters,[1 numClusters]); % Final state order
Xo_ind = repelem(1:numClusters,numClusters); % Paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % Isolate off diagonal from linearized transition probabilities

% Initialize values for the loop
x0 = Averaged_CentroidsZS_TS(:,Xo_ind);
xf = Averaged_CentroidsZS_TS(:,Xf_ind);
E_Full_T = NaN(nT,numClusters^2);

for i = 1:nT
    T=T_range(i);
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    ndims(Anorm);
    E_Full_T(i,:) = MIN_CONTROL_ENERGY_new(Anorm, WcI, x0, xf, T, false); % compute minimum control energy for each state transition
    disp(E_Full_T(i:6,1:6));
    
    disp(['Current value of T: ', num2str(T)]);
    
end


E_Full_T_k6_ZS_TS = E_Full_T;
E_full_T_k6_SC_Diag0_c1 = E_Full_T;

% Save and load energy data
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilE.Full.T.k6.ZS.TS.mat','E_Full_T_k6_ZS_TS');
% E_full_grpavg_T_k6_zscoredTS = load('/Users/daraneumann/Downloads/Research/NCT/Results/.Full.T.k6.ZS.TS.mat').E_Full_T_k6_ZS_TS;
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilE.Full.T.k6.SC.Diag0.c1.mat','E_full_T_k6_SC_Diag0_c1');
% E_full_grpavg_T_k6_zscoredTS = load('/Users/daraneumann/Downloads/Research/NCT/Results/E.Full.T.k6.SC.Diag0.c1.mat').E_full_T_k6_SC_Diag0_c1;


%% CONTROL ENERGY FOR T = 0.501


% Locate the MATLAB function
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/control'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/dynamicsfxns'

% Initialize values
Centroids_PerPt_ZS_TS=load('/Users/daraneumann/Downloads/Research/NCT/Results/AprilCentroids.PerPt.ZS.TS.mat').Centroids_PerPt_ZS_TS;
numClusters = 6;
T = 0.501;
c = 1;
SC = load('/Users/daraneumann/Downloads/Research/NCT/Data/SC.Matrix.diag0.mat').labpcexport;
n = 554;

% Visualize SC map
SCdataforNCT1 = SC;
size(SCdataforNCT1)
imagesc(SCdataforNCT1(:,:,1))

% Run a loop for SC data for all brain regions
SCdataforNCT1 = zeros(n,200,200);

for s = 1:n
    SCdataforNCT1(s,:,:) = SC;
end
size(SCdataforNCT1)
imagesc(SCdataforNCT1(:,:,1))


% Compute minimum control energy
minimumcontrolenergy = zeros(n,36);

E_Full = NaN(n,numClusters^2);size(E_Full) % State-wise transition energy
E_Regional = NaN(n,numClusters^2,200);size(E_Regional) % Regional transition energy

for i = 1:n
    SCdataforNCT2 = SCdataforNCT1(i,:,:);
    A = reshape(SCdataforNCT2,200,200);
    Anorm = NORMALIZE(A,c); % Normalize A by maximum eigenvalue - eye(N) to make marginally stable
    
    % Define x0, xf, initial state, and final state as cluster centroids for each
    Xf_ind = repmat(1:numClusters,[1 numClusters]); % Final state order
    Xo_ind = repelem(1:numClusters,numClusters); % Paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
    onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
    offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % Isolate off diagonal from linearized transition probabilities
    
    centroids = reshape(Centroids_PerPt_ZS_TS(i,:,:),200,6);
    
    x0 = centroids(:,Xo_ind);
    xf = centroids(:,Xf_ind); % Allows each column of x0 and xf to represent state transitions
    WcI = GRAMIAN_FAST(Anorm, T); % Compute gramian inverse for control horizon T
    
    % E_Full = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % Compute minimum control energy for each state transition
    [E_Full(i,:), E_Regional(i,:,:)] = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % Compute minimum CE for each state and regional transition
    % plot(E_Full)
    % minimumcontrolenergy(i,:) = E_Full;
    % minimumcontrolenergy(k,i,:) = E_Full;
   
end

% Compute global and regional transition energy
globalTE = mean(E_Full,2);size(globalTE)
regionalTE = mean(E_Regional,2);size(regionalTE)

% Visualize the minimum control energies
minimumcontrolenergy = E_Full;
size(minimumcontrolenergy)
imagesc(minimumcontrolenergy)
colMeans = mean(minimumcontrolenergy, 1);

% Generate matrix
CEmatrix = reshape(colMeans,6,6);
imagesc(CEmatrix);colormap('default');
xticklabels({'Limbic (-)','SomMotor (-)','Visual (-)','Visual (+)','SomMotor (+)','Limbic (+)'});
yticklabels({'Limbic (-)','SomMotor (-)','Visual (-)','Visual (+)','SomMotor (+)','Limbic (+)'});
colorbar;

% Order the matrix to make network comparisons on the graph more intuitive
CEmatrix_ordered = CEmatrix([6,1,5,2,4,3],[3,4,2,5,1,6]); % First argument used for the y-axis, which has 6 at the origin. Change based on R
imagesc(CEmatrix_ordered);colormap('default');
xticklabels({'Limbic (+)','Limbic (-)','SomMotor (+)','SomMotor (-)','Visual (+)','Visual (-)'});
yticklabels({'Visual (-)','Visual (+)','SomMotor (-)','SomMotor (+)','Limbic (-)','Limbic (+)'});
colorbar;

% Save state-wise, global, and regional minimum control energy
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilMinimum.Control.Energy.mat','minimumcontrolenergy');
% minimumcontrolenergy = load('/Users/daraneumann/Downloads/Research/NCT/Results/Minimum.Control.Energy.mat').minimumcontrolenergy
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilGlobal.Minimum.CE.mat','globalTE');
% globalTE = load('/Users/daraneumann/Downloads/Research/NCT/Results/AprilGlobal.Minimum.CE.mat').globalTE
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilRegional.Minimum.CE.mat','regionalTE');
% regionalTE = load('/Users/daraneumann/Downloads/Research/NCT/Results/AprilRegional.Minimum.CE.mat').regionalTE
%% CONTROL ENERGY FOR T = 1


% Locate the MATLAB function
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/control'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/dynamicsfxns'

% Initialize values
Centroids_PerPt_ZS_TS=load('/Users/daraneumann/Downloads/Research/NCT/Results/AprilCentroids.PerPt.ZS.TS.mat').Centroids_PerPt_ZS_TS;numClusters = 6;
T = 1;
c = 1;
SC = load('/Users/daraneumann/Downloads/Research/NCT/Data/SC.Matrix.diag0.mat').labpcexport;
n = 554;

% Visualize SC map
SCdataforNCT1 = SC;
size(SCdataforNCT1)
imagesc(SCdataforNCT1(:,:,1))

% Run a loop for SC data for all brain regions
SCdataforNCT1 = zeros(n,200,200);

for s = 1:n
    SCdataforNCT1(s,:,:) = SC;
end
size(SCdataforNCT1)
imagesc(SCdataforNCT1(:,:,1))


% Compute minimum control energy
minimumcontrolenergy = zeros(n,36);

E_Full = NaN(n,numClusters^2);size(E_Full) % State-wise transition energy
E_Regional = NaN(n,numClusters^2,200);size(E_Regional) % Regional transition energy

for i = 1:n
    SCdataforNCT2 = SCdataforNCT1(i,:,:);
    A = reshape(SCdataforNCT2,200,200);
    Anorm = NORMALIZE(A,c); % Normalize A by maximum eigenvalue - eye(N) to make marginally stable
    
    % Define x0, xf, initial state, and final state as cluster centroids for each
    Xf_ind = repmat(1:numClusters,[1 numClusters]); % Final state order
    Xo_ind = repelem(1:numClusters,numClusters); % Paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
    onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
    offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % Isolate off diagonal from linearized transition probabilities
    
    centroids = reshape(Centroids_PerPt_ZS_TS(i,:,:),200,6);
    
    x0 = centroids(:,Xo_ind);
    xf = centroids(:,Xf_ind); % Allows each column of x0 and xf to represent state transitions
    WcI = GRAMIAN_FAST(Anorm, T); % Compute gramian inverse for control horizon T
    
    % E_Full = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % Compute minimum control energy for each state transition
    [E_Full(i,:), E_Regional(i,:,:)] = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % Compute minimum CE for each state and regional transition
    % plot(E_Full)
    % minimumcontrolenergy(i,:) = E_Full;
    % minimumcontrolenergy(k,i,:) = E_Full;
   
end

% Compute global and regional transition energy
globalTE = mean(E_Full,2);size(globalTE)
regionalTE = mean(E_Regional,2);size(regionalTE)

% Visualize the minimum control energies
minimumcontrolenergy = E_Full;
size(minimumcontrolenergy)
imagesc(minimumcontrolenergy)
colMeans = mean(minimumcontrolenergy, 1);

% Generate matrix
CEmatrix = reshape(colMeans,6,6);
imagesc(CEmatrix);colormap('default');
xticklabels({'Limbic (-)','SomMotor (-)','Visual (-)','Visual (+)','SomMotor (+)','Limbic (+)'});
yticklabels({'Limbic (-)','SomMotor (-)','Visual (-)','Visual (+)','SomMotor (+)','Limbic (+)'});
colorbar;

% Order the matrix to make network comparisons on the graph more intuitive
CEmatrix_ordered = CEmatrix([6,1,5,2,4,3],[3,4,2,5,1,6]); % First argument used for the y-axis, which has 6 at the origin. Change based on R
imagesc(CEmatrix_ordered);colormap('default');
xticklabels({'Limbic (+)','Limbic (-)','SomMotor (+)','SomMotor (-)','Visual (+)','Visual (-)'});
yticklabels({'Visual (-)','Visual (+)','SomMotor (-)','SomMotor (+)','Limbic (-)','Limbic (+)'});
colorbar;

% Save state-wise, global, and regional minimum control energy
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilMinimum.Control.Energy.T1.mat','minimumcontrolenergy');
% minimumcontrolenergy = load('/Users/daraneumann/Downloads/Research/NCT/Results/Minimum.Control.Energy.T1.mat').minimumcontrolenergy
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilGlobal.Minimum.CE.T1.mat','globalTE');
% globalTE = load('/Users/daraneumann/Downloads/Research/NCT/Results/Global.Minimum.CE.T1.mat').globalTE
 save('/Users/daraneumann/Downloads/Research/NCT/Results/AprilRegional.Minimum.CE.T1.mat','regionalTE');
% regionalTE = load('/Users/daraneumann/Downloads/Research/NCT/Results/Regional.Minimum.CE.T1.mat').regionalTE

%% BRAIN PLOTS FOR STATE-WISE TRANSITION ENERGY CENTROIDS T=0.501

% Make sure to run "CONTROL ENERGY FOR T = 0.501" first

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; %#ok<NASGU>


results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';
best_number_of_clusters = 6;


% Locate MATLAB function
whichatlas = 'schaefer200x7'; %#ok<NASGU>
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
cmap = colormap(plasma);

% Set scale between 1 and -1
data_min = min(min(centroids));
data_max = max(max(centroids));
centroids1 = centroids/max(abs(data_min),data_max);

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'

% Create images
for i = 1:best_number_of_clusters
    close all;
    data = centroids1(:,i);
    
    img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap);
    
    figure('Position', [0 0 2000 2000]);
    imshow(img);
    set(gca,'colormap',cmap);
    c = colorbar('SouthOutside', 'fontsize', 20);
    clim([-1 1]);
    % c.Label.String = 'Colorbar Label';
    % caxis([data_min data_max]);
    
    % Save images
    saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.Cluster.', num2str(i), '.png'));
    pause(2)
end

%% BRAIN PLOTS FOR STATE-WISE TRANSITION ENERGY CENTROIDS T=1

% Make sure to run "CONTROL ENERGY FOR T = 1" first

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; %#ok<NASGU>


results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';
best_number_of_clusters = 6;


% Locate MATLAB function
whichatlas = 'schaefer200x7'; %#ok<NASGU>
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
cmap = colormap(plasma);

% Set scale between 1 and -1
data_min = min(min(centroids));
data_max = max(max(centroids));
centroids1 = centroids/max(abs(data_min),data_max);

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'

% Create images
for i = 1:best_number_of_clusters
    close all;
    data = centroids1(:,i);
    
    img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap);
    
    figure('Position', [0 0 2000 2000]);
    imshow(img);
    set(gca,'colormap',cmap);
    c = colorbar('SouthOutside', 'fontsize', 20);
    clim([-1 1]);
    % c.Label.String = 'Colorbar Label';
    % caxis([data_min data_max]);
    
    % Save images
    saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.Cluster.T1.', num2str(i), '.png'));
    pause(2)
end

%% (OLD) BRAINPLOTS FOR REGIONAL TRANSITION ENERGY CENTROIDS

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; 

estimates_status = load('/Users/daraneumann/Downloads/Research/NCT/Data/CE.Estimates.Status.mat').labpcexport;
estimates_sex = load('/Users/daraneumann/Downloads/Research/NCT/Data/CE.Estimates.Sex.mat').labpcexport;
results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Locate MATLAB function
whichatlas = 'schaefer200x7'; 
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
cmap = colormap(plasma);

% Set scale between 1 and -1
data_min_status = min(min(estimates_status));
data_max_status = max(max(estimates_status));
centroids_status = estimates_status/max(abs(data_min_status),data_max_status);

data_min_sex = min(min(estimates_sex));
data_max_sex = max(max(estimates_sex));
centroids_sex = estimates_sex/max(abs(data_min_sex),data_max_sex);

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'

% Create images for status
close all;
data = centroids_status;

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap);

figure('Position', [0 0 2000 2000]);
imshow(img);
c = colorbar('SouthOutside', 'fontsize', 20);
c.Label.String = 'Colorbar Label';
set(gca,'colormap',cmap);
caxis([-1 1]); %#ok<CAXIS>
% caxis([data_min data_max]);
    
% Add labels
% title('TT', 'fontsize',18);
% annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',20,'color','white')
% annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',20,'color','white')
% annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center','backgroundcolor','black')
% annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')
% set(gcf, 'Position', [0 0 2000 2000]);

% Save image for status
saveas(gcf, strcat(results_dir, 'Centroids.Status.', num2str(i), '.png'));
pause(2)


% Create images for sex
close all;
data = centroids_sex;

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap);

figure('Position', [0 0 2000 2000]);
imshow(img);
c = colorbar('SouthOutside', 'fontsize', 20);
c.Label.String = 'Colorbar Label';
set(gca,'colormap',cmap);
caxis([-1 1]); %#ok<CAXIS>
% caxis([data_min data_max]);
    
% Add labels
% title('TT', 'fontsize',18);
% annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',20,'color','white')
% annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',20,'color','white')
% annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center','backgroundcolor','black')
% annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',20,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')
% set(gcf, 'Position', [0 0 2000 2000]);

% Save image for sex
saveas(gcf, strcat(results_dir, 'Centroids.Sex.', num2str(i), '.png'));
pause(2)


%% BRAIN PLOTS FOR REGIONAL ENTROPIES BY STATUS

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs;


% Locate MATLAB function
whichatlas = 'schaefer200x7';
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
cmap = colormap(viridis);
addpath '/Users/daraneumann/Downloads/Research/NCT/Code'
RdBu = csvread('RdBu_colormap.csv', 1, 0);

% Load data
estimates_entropy = load('/Users/daraneumann/Downloads/Research/NCT/Data/Entropy.Estimates.Status.mat').labpcexport;
estimates_entropy_t = load('/Users/daraneumann/Downloads/Research/NCT/Data/Entropy.T.mat').labpcexport;
average_entropy_HC = load('/Users/daraneumann/Downloads/Research/NCT/Data/Average.Regional.Entropy.HC.mat').labpcexport;
average_entropy_MCI = load('/Users/daraneumann/Downloads/Research/NCT/Data/Average.Regional.Entropy.MCI.mat').labpcexport;
results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Set scale for the plots showing averages (must be the same scale)
entropy_min = min(min(average_entropy_HC),min(average_entropy_MCI)); % Find lowest of the minimums to set as the lower bound of the graph
entropy_max = max(max(average_entropy_HC),max(average_entropy_MCI)); % Find highest of the maximums to set as the upper bound of the graph

% Create images for HC
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'
data = average_entropy_HC;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[entropy_min 0.75]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[entropy_min 0.75]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([entropy_min entropy_max]); % Adjust based on entropy_min and entropy_max
% c.Label.String = 'Colorbar Label';

% Save image for HC
saveas(gcf, strcat(results_dir, 'Brain.Plot.Entropy.HC.png'));
pause(2)


% Create images for MCI
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'
data = average_entropy_MCI;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[entropy_min 0.75]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[entropy_min 0.75]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([entropy_min entropy_max]) % Adjust based on entropy_min and entropy_max
% c.Label.String = 'Colorbar Label';

% Save image for MCI
saveas(gcf, strcat(results_dir, 'Brain.Plot.Entropy.MCI.png'));
pause(2)


% Create images for t-statistic comparisons
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'
data = estimates_entropy;
cmap = RdBu;
%cmap = colormap(viridis); % Select color map

entropy_t_max = max(estimates_entropy);
entropy_t_min = min(estimates_entropy);

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[-0.1693 0.1693]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[-0.1693 0.1693]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';

% Save image for t-statistic comparisons
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.Entropy.Comparison.png'));
pause(2)

%% BRAIN PLOTS FOR TRANSITION ENERGIES BY STATUS T=0.501

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; 

% Import data
estimates_CE = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilCE.Estimates.Status.mat').labpcexport;
average_CE_HC = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilAverage.Regional.CE.HC.mat').labpcexport;
average_CE_MCI = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilAverage.Regional.CE.MCI.mat').labpcexport;
results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Locate MATLAB function
whichatlas = 'schaefer200x7'; 
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/slanCM/slanCM'
load('slanCM_Data.mat')
addpath '/Users/daraneumann/Downloads/Research/NCT/Code'
RdBu = csvread('RdBu_colormap.csv', 1, 0);

% Set scale for the plots showing averages (must be the same scale)
CE_min = min(min(average_CE_HC),min(average_CE_MCI)); % Find lowest of the minimums to set as the lower bound of the graph
CE_max = max(max(average_CE_HC),max(average_CE_MCI)); % Find highest of the maximums to set as the upper bound of the graph

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'

% Create images for HC
close all;
data = average_CE_HC;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0 CE_max]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0 CE_max]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on CE_min and CE_max
% c.Label.String = 'Colorbar Label';
    
% Save image for HC
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.CE.HC.png'));
pause(2)


% Create images for MCI


close all;
data = average_CE_MCI;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0 CE_max]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0 CE_max]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on CE_min and CE_max
% c.Label.String = 'Colorbar Label';
    
% Save image for MCI
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.CE.MCI.png'));
pause(2)


% Create images for t-statistic comparisons
close all;
data = estimates_CE;
cmap = colormap(turbo); % Select color map
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};

CE_t_max = max(estimates_CE);
CE_t_min = min(estimates_CE);

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-2.646 2.646]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-2.646 2.646]); % Can change clim based on CE_t_max and CE_t_min
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for t-statistic comparisons
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.CE.Comparison.png'));
pause(2)

%% BRAIN PLOTS FOR TRANSITION ENERGIES BY STATUS T=1

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; 

% Import data
estimates_CE_T1 = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilCE.Estimates.Status.T1.mat').labpcexport;
average_CE_HC_T1 = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilAverage.Regional.CE.HC.T1.mat').labpcexport;
average_CE_MCI_T1 = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilAverage.Regional.CE.MCI.T1.mat').labpcexport;
results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Locate MATLAB function
whichatlas = 'schaefer200x7'; 
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/slanCM/slanCM'
load('slanCM_Data.mat')
addpath '/Users/daraneumann/Downloads/Research/NCT/Code'
RdBu = csvread('RdBu_colormap.csv', 1, 0);

% Set scale for the plots showing averages (must be the same scale)
CE_min = min(min(average_CE_HC_T1),min(average_CE_MCI_T1)); % Find lowest of the minimums to set as the lower bound of the graph
CE_max = max(max(average_CE_HC_T1),max(average_CE_MCI_T1)); % Find highest of the maximums to set as the upper bound of the graph

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'

% Create images for HC
close all;
data = average_CE_HC_T1;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0 CE_max]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0 CE_max]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on CE_min and CE_max
% c.Label.String = 'Colorbar Label';
    
% Save image for HC
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.CE.HC.T1.png'));
pause(2)


% Create images for MCI


close all;
data = average_CE_MCI_T1;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0 CE_max]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0 CE_max]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on CE_min and CE_max
% c.Label.String = 'Colorbar Label';
    
% Save image for MCI
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.CE.MCI.T1.png'));
pause(2)


% Create images for t-statistic comparisons
close all;
data = estimates_CE_T1;
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};
cmap = colormap(turbo); % Select color map

CE_t_max = max(estimates_CE_T1);
CE_t_min = min(estimates_CE_T1);

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-1.6346 1.6346]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-1.6346 1.6346]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for t-statistic comparisons
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.CE.Comparison.T1.png'));
pause(2)

%% BRAIN PLOTS FOR PLAQUE BY STATUS

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; 

% Import data
estimates_plaque = load('/Users/daraneumann/Downloads/Research/NCT/Data/Plaque.Estimates.Status.mat').labpcexport;
estimates_plaque_t = load('/Users/daraneumann/Downloads/Research/NCT/Data/t.Estimates.Plaque.mat').labpcexport;
average_plaque_HC = load('/Users/daraneumann/Downloads/Research/NCT/Data/Average.Regional.Plaque.HC.mat').labpcexport;
average_plaque_MCI = load('/Users/daraneumann/Downloads/Research/NCT/Data/Average.Regional.Plaque.MCI.mat').labpcexport;
results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Locate MATLAB function
whichatlas = 'schaefer200x7'; 
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/slanCM/slanCM'
load('slanCM_Data.mat')
addpath '/Users/daraneumann/Downloads/Research/NCT/Code'
RdBu = csvread('RdBu_colormap.csv', 1, 0);

% Set scale for the plots showing averages (must be the same scale)
plaque_min = min(min(average_plaque_HC),min(average_plaque_MCI)); % Find lowest of the minimums to set as the lower bound of the graph
plaque_max = max(max(average_plaque_HC),max(average_plaque_MCI)); % Find highest of the maximums to set as the upper bound of the graph

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'

% Create images for HC
close all;
data = average_plaque_HC;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0.5 2]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0.5 2]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on plaque_min and plaque_max
% c.Label.String = 'Colorbar Label';
    
% Save image for HC
saveas(gcf, strcat(results_dir, 'Brain.Plot.Plaque.HC.png'));
pause(2)


% Create images for MCI
close all;
data = average_plaque_MCI;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0.5 2]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0.5 2]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on plaque_min and plaque_max
% c.Label.String = 'Colorbar Label';
    
% Save image for MCI
saveas(gcf, strcat(results_dir, 'Brain.Plot.Plaque.MCI.png'));
pause(2)


% Create images for t-statistic comparisons
close all;
data = estimates_plaque;
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};
cmap = colormap(turbo); % Select color map

plaque_lm_max = max(estimates_plaque);
plaque_lm_min = min(estimates_plaque);

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-0.4573 0.4573]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-0.4573 0.4573]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for t-statistic comparisons
saveas(gcf, strcat(results_dir, 'Brain.Plot.Plaque.Comparison.png'));
pause(2)

%% BRAIN PLOTS FOR TAU BY STATUS

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; 

% Import data
estimates_tau = load('/Users/daraneumann/Downloads/Research/NCT/Data/Tau.Estimates.Status.mat').labpcexport;
estimates_tau_t = load('/Users/daraneumann/Downloads/Research/NCT/Data/t.Estimates.Tau.mat').labpcexport;
average_tau_HC = load('/Users/daraneumann/Downloads/Research/NCT/Data/Average.Regional.Tau.HC.mat').labpcexport;
average_tau_MCI = load('/Users/daraneumann/Downloads/Research/NCT/Data/Average.Regional.Tau.MCI.mat').labpcexport;
results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Locate MATLAB function
whichatlas = 'schaefer200x7'; 
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code'
RdBu = csvread('RdBu_colormap.csv', 1, 0);

% Set scale for the plots showing averages (must be the same scale)
tau_min = min(min(average_tau_HC),min(average_tau_MCI)); % Find lowest of the minimums to set as the lower bound of the graph
tau_max = max(max(average_tau_HC),max(average_tau_MCI)); % Find highest of the maximums to set as the upper bound of the graph

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'

% Create images for HC
close all;
data = average_tau_HC;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0.75 2.75]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0.75 2.75]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on plaque_min and plaque_max
% c.Label.String = 'Colorbar Label';
    
% Save image for HC
saveas(gcf, strcat(results_dir, 'Brain.Plot.Tau.HC.png'));
pause(2)


% Create images for MCI
close all;
data = average_tau_MCI;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0.75 2.75]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[0.75 2.75]); % Can change clim, but must match above
c = colorbar('SouthOutside', 'fontsize', 20);
% caxis([0 4000]) % Adjust based on plaque_min and plaque_max
% c.Label.String = 'Colorbar Label';
    
% Save image for MCI
saveas(gcf, strcat(results_dir, 'Brain.Plot.Tau.MCI.png'));
pause(2)


% Create images for t-statistic comparisons
close all;
data = estimates_tau;
cmap = RdBu;
%cmap = colormap(viridis); % Select color map

tau_lm_max = max(estimates_tau);
tau_lm_min = min(estimates_tau);

%img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[0.0635 1.485]); % Can change clim
img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[-1.485 1.485]);

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[-1.485 1.485]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for t-statistic comparisons
saveas(gcf, strcat(results_dir, 'Brain.Plot.Tau.Comparison.png'));
pause(2)

%% BRAIN PLOT FOR TE VS ENTROPY

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; 

% Import data
cor_estimates_TE_entropy = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilTE.Entropy.Cor.Estimates.mat').labpcexport;
lm_estimates_TE_entropy = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilTE.Entropy.Estimates.mat').labpcexport;

results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Locate MATLAB function
whichatlas = 'schaefer200x7'; 
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/slanCM/slanCM' % Has the 'seismic' diverging color map
load('slanCM_Data.mat')
addpath '/Users/daraneumann/Downloads/Research/NCT/Code'
RdBu = csvread('RdBu_colormap.csv', 1, 0);

% Set scale for the plots (must be the same scale)
TE_entropy_min = min(lm_estimates_TE_entropy); % Find minimum to set as the lower bound of the graph
TE_entropy_max = max(lm_estimates_TE_entropy); % Find hthe maximum to set as the upper bound of the graph


addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'


% Create images for lm comparisons (TE and entropy)
close all;
data = lm_estimates_TE_entropy;
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};
cmap = colormap(turbo); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-0.475 0.475]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-4.75 4.75]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for lm comparisons
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.TE.Entropy.lm.png'));
pause(2)

%% BRAIN PLOTS FOR PET CORRELATIONS

% Locate MATLAB function
atlasblobs_list = load('/Users/daraneumann/Downloads/Research/NCT/Code/Other/AtlasBlobs.mat');
atlasblobs_list = atlasblobs_list.atlasblobs; 

% Import data
estimates_TE_plaque = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilTE.Plaque.Cor.Estimates.mat').labpcexport;
estimates_TE_tau = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilTE.Tau.Cor.Estimates.mat').labpcexport;
estimates_entropy_plaque = load('/Users/daraneumann/Downloads/Research/NCT/Data/Entropy.Plaque.Cor.Estimates.mat').labpcexport;
estimates_entropy_tau = load('/Users/daraneumann/Downloads/Research/NCT/Data/Entropy.Tau.Cor.Estimates.mat').labpcexport;

results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

% Locate MATLAB function
whichatlas = 'schaefer200x7'; 
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/ejc_bs_code/plottingfxns/Colormaps/Colormaps (5)/Colormaps'
addpath '/Users/daraneumann/Downloads/Research/NCT/Code/slanCM/slanCM' % Has the 'seismic' diverging color map
load('slanCM_Data.mat')
addpath '/Users/daraneumann/Downloads/Research/NCT/Code'
RdBu = csvread('RdBu_colormap.csv', 1, 0);

% Set scale for the plots (must be the same scale)
TE_plaque_min = min(estimates_TE_plaque); % Find minimum to set as the lower bound of the graph
TE_plaque_max = max(estimates_TE_plaque); % Find hthe maximum to set as the upper bound of the graph

TE_tau_min = min(estimates_TE_tau);
TE_tau_max = max(estimates_TE_tau);

entropy_plaque_min = min(estimates_entropy_plaque);
entropy_plaque_max = max(estimates_entropy_plaque);

entropy_tau_min = min(estimates_entropy_tau);
entropy_tau_max = max(estimates_entropy_tau);

%overall_min = min(min(estimates_TE_plaque),min(estimates_TE_tau),min(estimates_entropy_plaque),min(estimates_entropy_tau));
%overall_max = max(max(estimates_TE_plaque),max(estimates_TE_tau),max(estimates_entropy_plaque),max(estimates_entropy_tau));

addpath '/Users/daraneumann/Downloads/Research/NCT/Code/AtlasBlobs'


% Create images for comparisons (TE and plaque)
close all;
data = estimates_TE_plaque;
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};
cmap = colormap(turbo); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for comparisons
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.TE.Plaque.Spearman.png'));
pause(2)



% Create images for comparisons (TE and tau)
close all;
data = estimates_TE_tau;
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};
cmap = colormap(turbo); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for comparisons
saveas(gcf, strcat(results_dir, 'AprilBrain.Plot.TE.Tau.Spearman.png'));
pause(2)



% Create images for comparisons (entropy and plaque)
close all;
data = estimates_entropy_plaque;
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};
cmap = colormap(turbo); % Select color map


img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for comparisons
saveas(gcf, strcat(results_dir, 'Brain.Plot.Entropy.Plaque.Spearman.png'));
pause(2)



% Create images for comparisons (entropy and tau)
close all;
data = estimates_entropy_tau;
colmap = RdBu;
%colmap = slandarerCM(6).Colors{strcmpi(slandarerCM(6).Names, 'seismic')};
cmap = colormap(turbo); % Select color map


img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[-0.4325 0.4325]); % Can change clim, but must match above
c = colorbar('WestOutside', 'fontsize', 20);
% c.Label.String = 'Colorbar Label';
    
% Save image for comparisons
saveas(gcf, strcat(results_dir, 'Brain.Plot.Entropy.Tau.Spearman.png'));
pause(2)

%% ADJUSTED MUTUAL INFORMATION

newTS = load('/Users/daraneumann/Downloads/Research/NCT/Data/AprilTimeSeriesZScored.mat').labpcexport;
results_dir = '/Users/daraneumann/Downloads/Research/NCT/Results/Figures/';

distanceMethod = 'correlation'; % Distance metric for clustering, we used correlation
nreps = 50;	% How many times to repeat clustering; will choose lowest error solution
maxI = 1000; % How many times to let kmeans try to converge before aborting rep
split = 30;

concTS = newTS; % The time series matrix you used to cluster at the very beginning of the analysis

[T,nparc] = size(concTS);


% Generate N partitions that will be compared pair-wise for mutual information

for numClusters = [5 6] % If undecided, run loop; if confident about k, you can run over only 1 value for numClusters
    disp(['Starting clusters number: ',num2str(numClusters)]);
    N = 10; % Number of partitions to create and compare
    parts = NaN(T,N); % Partitions will be stored here
    D = NaN(N,T,numClusters); % Distance matrices will be stored here
end
    
    %load centroids to initialize kmeans on (optional) 
    % MUST ADD <'Start', init> to the kmeans inputs to use this
%     load(['Partition_bp12_k',num2str(numClusters),'.mat'],'centroids'); %make sure file matches centroids you want to initialize from
%     init = NaN(numClusters,nparc,nreps);
%     for i=1:nreps
%         init(:,:,i) = centroids.';
%     end
    
    for i=1:N
        disp(['Clusters: ',num2str(numClusters),'. Starting kmeans number: ',num2str(i)]);
        [parts(:,i),~,~,D(i,:,:)] = kmeans(concTS,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    end
    
    % Calculate adjusted mutual information for every pair of partitions
    addpath '/Users/daraneumann/Downloads/Research/NCT/Code/EnergyLandscape/misc_code'
    ami_results = NaN(N,N);
    
    for i=1:N
        for j=1:N
            ami_results(i,j) = ami(parts(:,i),parts(:,j));
        end
    end
    
    % Assess
    [m,ind] = max(sum(ami_results,1)); % ind corresponds to the partition which has the highest mutual information with all other partitions
    partition = parts(:,ind); % Take partition that has most agreement with all other for further analysis
    
    % Plot
    f = figure;
    
    imagesc(ami_results); title(['Adjusted Mutual Information between Partitions k = ',num2str(numClusters)]); colorbar;
    axis square; set(gca,'FontSize',8);
    f.PaperUnits = 'inches';
    f.PaperSize = [6 6];
    f.PaperPosition = [0 0 6 6];
    saveas(f,fullfile(results_dir,['AprilAMI_bp',num2str(split),num2str(numClusters),'.pdf']));
    


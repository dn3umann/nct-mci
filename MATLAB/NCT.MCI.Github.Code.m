
%% ---- 1) CONFIGURATION ---- %%

% Set up relative file paths (see README for file structure)
base_dir = '/Users/daraneumann/Downloads/Research/NCT';
data_dir = fullfile(base_dir, 'Data');
matlab_dir = fullfile(base_dir, 'MATLAB');
packages_dir = fullfile(matlab_dir, 'Packages');
results_dir = fullfile(base_dir, 'Results');
figures_dir = fullfile(results_dir, 'Figures');





%% ---- 2) COMPUTE BEST NUMBER OF CLUSTERS ---- %%

% Load the z-scored TS data from R
TS = load(fullfile(matlab_dir,"TimeSeriesZScored.mat")).labpcexport;

% Set parameters
distanceMethod = 'correlation'; % distance metric
nreps = 50; % number of times to repeat clustering; lowest error solution will be chosen
maxI = 1000; % tolerance for k-means convergence attempts
split = 10; % number of splits


% Compute ideal number of clusters between 2 and 10
for numClusters = 2:10

    disp(['Starting K = ',num2str(numClusters),' Split = ',num2str(split)]); % progress update for when a cluster is being started

    [partition,~,sumd] = kmeans(TS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    save(fullfile(results_dir,"Clusters",['KMeans',num2str(split),'.K',num2str(numClusters),'.mat']),'partition','sumd') % save data for each cluster
    
    disp(['Finished K = ',num2str(numClusters),' Split = ',num2str(split)]); % progress update for when a cluster is finished
    
    clear partition
    clear sumd 
end

disp('Done computing clusters.'); % progress update for end of loop





%% ---- 3) VARIANCE EXPLAINED BY CLUSTERING ---- %%

% Locate MATLAB packages
addpath(packages_dir);

% Initialize values for the loop
k_range = 2:10; % number of clusters used
VarianceExplained = zeros(length(k_range),1);

% Run the loop to calculate centroids and variance explained for each cluster number
for numClusters = k_range
    disp(['K = ',num2str(numClusters)])

    load(fullfile(results_dir,"Clusters",['KMeans',num2str(split),'.K',num2str(numClusters),'.mat']); % import data from last section
   
    kClusterCentroids = GET_CENTROIDS(TS,partition,numClusters);
    VarianceExplained(numClusters - 1) = VAREXPLAINED(TS,partition,kClusterCentroids,numClusters);
end

% Create figure for variance explained
f=figure;
plot(k_range,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
set(f, 'PaperUnits', 'inches');
x_width=6 ;y_width=6;
set(f, 'PaperPosition', [0 0 x_width y_width]);
saveas(f,fullfile(figures_dir,"Variance.Explained.Graph.png"));


% Create figure for gain in variance explained
f = figure;
plot(k_range(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
set(f, 'PaperUnits', 'inches');
x_width=6 ;y_width=6;
set(f, 'PaperPosition', [0 0 x_width y_width]);
saveas(f,fullfile(figures_dir,"Gain.Variance.Explained.Graph.png"));






%% ---- 4) PARTITIONS AND CENTROIDS ---- %%

% Set parameters based of previous results
numClusters = 6; % enter best number of clusters
distanceMethod = 'correlation'; % distance metric
nreps = 50; % number of times to repeat clustering (should be the same as previous)

% Compute partitions for the z-scored time series
[partitions,~,sumd] = kmeans(TS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',1000);
size(partitions)

% Compute centroids for the partition
centroids = GET_CENTROIDS(TS,partitions,numClusters);

% Save partition and centroid data for chosen k
save(fullfile(matlab_dir,['Partitions.K',num2str(numClusters),'.mat']),"partitions")
save(fullfile(matlab_dir,'Centroids.k.mat'),"centroids")





%% ---- 5) ADJUSTED MUTUAL INFORMATION ---- %%

distanceMethod = 'correlation'; % distance metric for clustering
nreps = 50;	% number of times to repeat clustering; will choose lowest error solution
maxI = 1000; % number of times to let kmeans try to converge before aborting rep
split = 30;
[T,nparc] = size(TS);


% Generate N partitions that will be compared pair-wise for mutual information
for numClusters = [5 6] % if undecided, run loop; if confident about k, you can run over only 1 value for numClusters
    disp(['Starting clusters number: ',num2str(numClusters)]);
    N = 10; % number of partitions to create and compare
    parts = NaN(T,N); % partitions will be stored here
    D = NaN(N,T,numClusters); % distance matrices will be stored here
end


for i=1:N
    disp(['Clusters: ',num2str(numClusters),'. Starting kmeans number: ',num2str(i)]);
    [parts(:,i),~,~,D(i,:,:)] = kmeans(TS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);
end
    

% Calculate adjusted mutual information for every pair of partitions

ami_results = NaN(N,N);
    
for i=1:N
    for j=1:N
        ami_results(i,j) = ami(parts(:,i),parts(:,j));
    end
end
    
% Assess the results
[m,ind] = max(sum(ami_results,1)); % ind corresponds to the partition which has the highest mutual information with all other partitions
partition = parts(:,ind); % take partition that has most agreement with all other for further analysis

% Plot and save the figure
f = figure;
imagesc(ami_results); title(['Adjusted Mutual Information between Partitions k = ',num2str(numClusters)]); colorbar;
axis square; set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [6 6];
f.PaperPosition = [0 0 6 6];
saveas(f,fullfile(figures_dir,['AMI_bp',num2str(split),num2str(numClusters),'.pdf']));
    




%% ---- 6) TIME SERIES LENGTH ---- %%

% Import time series length data from R and load in previous partition data
length_timeseries = load(fullfile(matlab_dir,'LengthTS.mat')).labpcexport;
load(fullfile(matlab_dir,['Partitions.K',num2str(numClusters),'.mat']));

% Identify best number of clusters
best_number_of_clusters = max(partitions);
disp(best_number_of_clusters)

% Assign subject IDs based on time series length
n = 554;
subjID = [];
for i = 1:n
    leng = length_timeseries(i);
    subjID = [subjID; ones(leng, 1)*i]; %#ok<AGROW>
end





%% ---- 7) TRANSITION PROBABILITIES ---- %%

% Initialize values
best_number_of_clusters = 6;
n = 554;
[transProbs,transitionProbabilityMats,numTransitions] = GET_TRANS_PROBS(partitions,subjID,best_number_of_clusters);
tran_prob = cell(n,1);

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
save(fullfile(results_dir,"Transition.Prob.mat","tran_prob"));





%% ---- 8) CENTROIDS FOR ALL SUBJECTS ---- %%

% Initialize values
numClusters = 6; % adjust based on previous result
n = 554;

% Run a loop to compile centroids for all subjects
centroids_per_subject = zeros(n,200,numClusters);

for subjectno = 1:n
    row = find(subjID==subjectno);
    centroids_per_subject(subjectno,:,:) = GET_CENTROIDS(TS(row,:),partitions(row),numClusters);
end


% Validate size of the matrix
size(centroids_per_subject)

% Confirm that the map is complete
imagesc(reshape(centroids_per_subject(1,:,:),200,numClusters));

% Create a map of the averaged centroids across all subjects
averaged_centroids = reshape(mean(centroids_per_subject(1:n,:,:)),200,numClusters);imagesc(averaged_centroids);

% Save values of centroids for all subjects
save(fullfile(results_dir,"Centroids.Per.Subject.mat","centroids_per_subject"));





%% ---- 9) FIND THE BEST T ---- %%

% Load SC data from R
SC = load(fullfile(matlab_dir,"SC.Matrix.diag0.mat")).labpcexport;

% Initialize values
nperms = 1000;
numClusters = 6; % adjust based on previous results

% Set parameters
c = 1; % time scale parameter based on values from paper
Anorm = NORMALIZE(SC,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable
T_range = 0.001:0.5:10; nT = length(T_range); % choose range of T values to test
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

% Set up variables and run the loop
x0 = averaged_centroids(:,Xo_ind);
xf = averaged_centroids(:,Xf_ind);
E_Full_T = NaN(nT,numClusters^2);

for i = 1:nT
    T=T_range(i);
    WcI = GRAMIAN_FAST(Anorm,T); % compute gramian inverse for control horizon T
    ndims(Anorm);
    E_Full_T(i,:) = MIN_CONTROL_ENERGY_new(Anorm,WcI,x0,xf,T,false); % compute minimum control energy for each state transition
    disp(E_Full_T(i:6,1:6));
    
    disp(['Current value of T: ',num2str(T)]);
    
end

% Save energy data
save(fullfile(results_dir,"E.Full.T.mat","E_Full_T"));





%% ---- 10) TRANSITION ENERGY ---- %%

% Load data and initialize values
SC = load(fullfile(matlab_dir,"SC.Matrix.diag0.mat")).labpcexport;
centroids_per_subject = load(fullfile(results_dir,"Centroids.Per.Subject.mat")).labpcexport;
numClusters = 6; % adjust based on previous results
T = 0.501; % adjust based on previous results
c = 1;
n = 554;

% Visualize SC map
imagesc(SC(:,:,1))

% Get SC data for all brain regions
SC_initial = zeros(n,200,200);

for s = 1:n
    SC_initial(s,:,:) = SC;
end

size(SC_initial)
imagesc(SC_initial(:,:,1))

% Compute minimum transition energy
minimumTE = zeros(n,36); % second argument in zeros = k^2
E_Full = NaN(n,numClusters^2);size(E_Full) % variable for global transition energy
E_Regional = NaN(n,numClusters^2,200);size(E_Regional) % variable for regional transition energy

for i = 1:n
    SC_new = SC_initial(i,:,:);
    A = reshape(SC_new,200,200); % third and fourth arguments should be equal to number of regions
    Anorm = NORMALIZE(A,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable
    
    % Define x0, xf, initial state, and final state as cluster centroids for each
    Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
    Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
    onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
    offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonals from linearized transition probabilities
    
    centroids_TE = reshape(centroids_per_subject(i,:,:),200,6);
    
    x0 = centroids_TE(:,Xo_ind);
    xf = centroids_TE(:,Xf_ind); % allows each column of x0 and xf to represent state transitions
    WcI = GRAMIAN_FAST(Anorm,T); % compute gramian inverse for control horizon T
    
    % Compute minimum TE for state and regional transitions
    [E_Full(i,:), E_Regional(i,:,:)] = MIN_CONTROL_ENERGY(Anorm,WcI,x0,xf,T,false);
   
end

% Compute global and regional transition energy
globalTE = mean(E_Full,2);size(globalTE)
regionalTE = mean(E_Regional,2);size(regionalTE)

% Visualize the minimum TE
minimumTE = E_Full;
imagesc(minimumTE)
colMeans = mean(minimumTE, 1);

% Generate matrix for the networks determined in R
TEmatrix = reshape(colMeans,6,6);
imagesc(CEmatrix);colormap('default');
xticklabels({'Limbic (-)','SomMotor (-)','Visual (-)','Visual (+)','SomMotor (+)','Limbic (+)'});
yticklabels({'Limbic (-)','SomMotor (-)','Visual (-)','Visual (+)','SomMotor (+)','Limbic (+)'});
colorbar;

% Order the matrix to make network comparisons on the graph more intuitive
TEmatrix_ordered = TEmatrix([6,1,5,2,4,3],[3,4,2,5,1,6]); % first argument used for the y-axis, which has 6 at the origin. Change based on R
imagesc(TEmatrix_ordered);colormap('default');
xticklabels({'Limbic (+)','Limbic (-)','SomMotor (+)','SomMotor (-)','Visual (+)','Visual (-)'});
yticklabels({'Visual (-)','Visual (+)','SomMotor (-)','SomMotor (+)','Limbic (-)','Limbic (+)'});
colorbar;

% Save state-wise, global, and regional minimum transition energy
save(fullfile(results_dir,"State.Minimum.TE.mat","minimumTE");
save(fullfile(results_dir,"Global.Minimum.TE.mat","globalTE");
save(fullfile(results_dir,"Regional.Minimum.TE.mat","regionalTE");





%% ---- 11) BRAIN PLOTS FOR CENTROIDS ---- %%

%Initialize data
best_number_of_clusters = 6;
cmap = colormap(plasma); % choose colormap for the brain plots

% Set scale between 1 and -1
data_min = min(min(centroids_TE));
data_max = max(max(centroids_TE));
centroids1 = centroids_TE/max(abs(data_min),data_max);

% Create images for each cluster
for i = 1:best_number_of_clusters
    close all;
    data = centroids1(:,i);
    
    img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap);
    
    figure('Position', [0 0 2000 2000]);
    imshow(img);
    set(gca,'colormap',cmap);
    c = colorbar('SouthOutside','fontsize',20); % color bar position and size
    clim([-1 1]);

    % Save images
    filename = fullfile(figures_dir,['Brain.Plot.Cluster.',num2str(i),'.png']);
    saveas(gcf,filename);
    
    pause(2)
end





%% ---- 12) BRAIN PLOTS FOR REGIONAL TRANSITION ENERGY ---- %%

% Import data from R
TE_estimates = load(fullfile(data_dir,"TE.Estimates.Status.mat")).labpcexport;
TE_HC = load(fullfile(data_dir,"Regional.TE.HC.mat")).labpcexport;
TE_MCI = load(fullfile(data_dir,"Regional.TE.HC.mat")).labpcexport;

% Import RdBu color map
RdBu = csvread(fullfile(packages_dir,'RdBu_colormap.csv'),1,0);

% Find minimums and maximums of HC, MCI, and comparisons to set scale
max_TE_status = max([TE_HC(:); TE_MCI(:)]); % find the maximum to set as the upper bound of the graph
min_TE_status = min([TE_HC(:); TE_MCI(:)]); % find minimum to set as the lower bound of the graph
max_TE_estimates = max(TE_estimates); % find the maximum to set as the upper bound of the graph
min_TE_estimates = min(TE_estimates); % find minimum to set as the lower bound of the graph



% Create brain plots for HC only
close all;
data = TE_HC;
cmap = colormap(plasma); % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_TE_status max_TE_status]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_TE_status max_TE_status]); % ensure clim here matches clim in line 403
c = colorbar('SouthOutside','fontsize', 20); % color bar settings
    
% Save image for HC
saveas(gcf,strcat(figures_dir,'Brain.Plot.TE.HC.png'));
pause(2)



% Create brain plots for MCI only
close all;
data = TE_MCI;
cmap = colormap(plasma); % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_TE_status max_TE_status]); % ensure clim matches HC plots

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_TE_status max_TE_status]); % ensure clim matches HC plots
c = colorbar('SouthOutside','fontsize',20); % color bar settings
    
% Save image for MCI
saveas(gcf,strcat(figures_dir,'Brain.Plot.TE.MCI.png'));
pause(2)



% Create brain plots for lm estimates
close all;
data = TE_estimates;
cmap = RdBu; % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_TE_estimates max_TE_estimates]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[min_TE_estimates max_TE_estimates]); % ensure clim matches clim in line 439
c = colorbar('WestOutside','fontsize',20); % color bar settings
    
% Save image for lm estimates
saveas(gcf,strcat(figures_dir,'Brain.Plot.TE.Comparison.png'));
pause(2)





%% ---- 13) BRAIN PLOTS FOR REGIONAL ENTROPY ---- %%

% Import data from R
entropy_estimates = load(fullfile(data_dir,"Entropy.Estimates.Status.mat")).labpcexport;
entropy_HC = load(fullfile(data_dir,"Regional.Entropy.HC.mat")).labpcexport;
entropy_MCI = load(fullfile(data_dir,"Regional.Entropy.HC.mat")).labpcexport;

% Find minimums and maximums of HC, MCI, and comparisons to set scale
max_entropy_status = max([entropy_HC(:); entropy_MCI(:)]); % find the maximum to set as the upper bound of the graph
min_entropy_status = min([entropy_HC(:); entropy_MCI(:)]); % find minimum to set as the lower bound of the graph
max_entropy_estimates = max(entropy_estimates); % find the maximum to set as the upper bound of the graph
min_entropy_estimates = min(entropy_estimates); % find minimum to set as the lower bound of the graph



% Create brain plots for HC only
close all;
data = entropy_HC;
cmap = colormap(plasma); % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_entropy_status max_entropy_status]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_entropy_status max_entropy_status]); % ensure clim here matches clim in line 472
c = colorbar('SouthOutside','fontsize',20); % color bar settings
    
% Save image for HC
saveas(gcf,strcat(figures_dir,'Brain.Plot.Entropy.HC.png'));
pause(2)



% Create brain plots for MCI only
close all;
data = entropy_MCI;
cmap = colormap(plasma); % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_entropy_status max_entropy_status]); % ensure clim matches HC plots

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_entropy_status max_entropy_status]); % ensure clim matches HC plots
c = colorbar('SouthOutside','fontsize',20); % color bar settings
    
% Save image for MCI
saveas(gcf,strcat(figures_dir,'Brain.Plot.Entropy.MCI.png'));
pause(2)



% Create brain plots for lm estimates
close all;
data = entropy_estimates;
cmap = RdBu; % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_entropy_estimates max_entropy_estimates]); % can change clim if needed

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[min_entropy_estimates max_entropy_estimates]); % ensure clim matches clim in line 508
c = colorbar('WestOutside', 'fontsize', 20); % Color bar settings
    
% Save image for lm estimates
saveas(gcf, strcat(figures_dir, 'Brain.Plot.Entropy.Comparison.png'));
pause(2)





%% ---- 14) BRAIN PLOTS FOR REGIONAL TAU ---- %%

% Import data from R
tau_estimates = load(fullfile(data_dir,"Tau.Estimates.Status.mat")).labpcexport;
tau_HC = load(fullfile(data_dir,"Regional.Tau.HC.mat")).labpcexport;
tau_MCI = load(fullfile(data_dir,"Regional.Tau.HC.mat")).labpcexport;

% Find minimums and maximums of HC, MCI, and comparisons to set scale
max_tau_status = max([tau_HC(:); tau_MCI(:)]); % find the maximum to set as the upper bound of the graph
min_tau_status = min([tau_HC(:); tau_MCI(:)]); % find minimum to set as the lower bound of the graph
max_tau_estimates = max(tau_estimates); % find the maximum to set as the upper bound of the graph
min_tau_estimates = min(tau_estimates); % find minimum to set as the lower bound of the graph



% Create brain plots for HC only
close all;
data = tau_HC;
cmap = colormap(plasma); % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_tau_status max_tau_status]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_tau_status max_tau_status]); % ensure clim here matches clim in line 541
c = colorbar('SouthOutside','fontsize',20); % color bar settings
    
% Save image for HC
saveas(gcf,strcat(figures_dir,'Brain.Plot.Tau.HC.png'));
pause(2)



% Create brain plots for MCI only
close all;
data = tau_MCI;
cmap = colormap(plasma); % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_tau_status max_tau_status]); % ensure clim matches HC plots

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_tau_status max_tau_status]); % ensure clim matches HC plots
c = colorbar('SouthOutside','fontsize',20); % color bar settings
    
% Save image for MCI
saveas(gcf,strcat(figures_dir,'Brain.Plot.Tau.MCI.png'));
pause(2)



% Create brain plots for lm estimates
close all;
data = tau_estimates;
cmap = RdBu; % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_tau_estimates max_tau_estimates]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[min_tau_estimates max_tau_estimates]); % ensure clim matches clim in line 577
c = colorbar('WestOutside','fontsize',20); % color bar settings
    
% Save image for lm estimates
saveas(gcf,strcat(figures_dir,'Brain.Plot.Tau.Comparison.png'));
pause(2)





%% ---- 15) BRAIN PLOTS FOR REGIONAL AMYLOID BETA PLAQUE ---- %%

% Import data from R
plaque_estimates = load(fullfile(data_dir,"Plaque.Estimates.Status.mat")).labpcexport;
plaque_HC = load(fullfile(data_dir,"Regional.Plaque.HC.mat")).labpcexport;
tau_MCI = load(fullfile(data_dir,"Regional.Plaque.HC.mat")).labpcexport;

% Find minimums and maximums of HC, MCI, and comparisons to set scale
max_plaque_status = max([plaque_HC(:); tau_MCI(:)]); % find the maximum to set as the upper bound of the graph
min_plaque_status = min([plaque_HC(:); tau_MCI(:)]); % find minimum to set as the lower bound of the graph
max_plaque_estimates = max(plaque_estimates); % find the maximum to set as the upper bound of the graph
min_plaque_estimates = min(plaque_estimates); % find minimum to set as the lower bound of the graph



% Create brain plots for HC only
close all;
data = plaque_HC;
cmap = colormap(plasma); % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_plaque_status max_plaque_status]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_plaque_status max_plaque_status]); % Ensure clim here matches clim in line 610
c = colorbar('SouthOutside','fontsize',20); % Color bar settings
    
% Save image for HC
saveas(gcf,strcat(figures_dir,'Brain.Plot.Plaque.HC.png'));
pause(2)



% Create brain plots for MCI only
close all;
data = tau_MCI;
cmap = colormap(plasma); % select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_plaque_status max_plaque_status]); % ensure clim matches HC plots

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',cmap,'clim',[min_plaque_status max_plaque_status]); % ensure clim matches HC plots
c = colorbar('SouthOutside','fontsize',20); % color bar settings
    
% Save image for MCI
saveas(gcf, strcat(figures_dir, 'Brain.Plot.Plaque.MCI.png'));
pause(2)



% Create brain plots for lm estimates
close all;
data = plaque_estimates;
cmap = RdBu; % Select color map

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[min_plaque_estimates max_plaque_estimates]); % Can change clim if needed

figure('Position', [0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[min_plaque_estimates max_plaque_estimates]); % Ensure clim matches clim in line 401
c = colorbar('WestOutside', 'fontsize', 20); % Color bar settings
    
% Save image for lm estimates
saveas(gcf,strcat(figures_dir,'Brain.Plot.Plaque.Comparison.png'));
pause(2)





%% ---- 16) BRAIN PLOTS FOR REGIONAL ASSOCIATION BETWEEN TE AND ENTROPY ---- %%

% Import data from R
lm_estimates_TE_entropy = load(fullfile(results_dir, "TE.Entropy.lm.Estimates.mat")).labpcexport;

% Set scale for the plots (must be the same scale)
TE_entropy_min = min(lm_estimates_TE_entropy); % find minimum to set as the lower bound of the graph
TE_entropy_max = max(lm_estimates_TE_entropy); % find the maximum to set as the upper bound of the graph

% Create brain plots for lm estimates
close all;
data = lm_estimates_TE_entropy;
cmap = RdBu;

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[TE_entropy_min TE_entropy_max]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[TE_entropy_min TE_entropy_max]); % ensure clim matches clim on line 673
c = colorbar('WestOutside','fontsize',20); % color bar settings

% Save brain plot image
saveas(gcf,strcat(figures_dir,'Brain.Plot.TE.Entropy.lm.png'));
pause(2)





%% ---- 17) BRAIN PLOTS FOR REGIONAL CORRELATIONS BETWEEN PET, TE, AND ENTROPY ---- %%

% Import data from R
cor_estimates_TE_plaque = load(fullfile(results_dir,"TE.Plaque.Cor.Estimates.mat")).labpcexport;
cor_estimates_TE_tau = load(fullfile(results_dir,"TE.Tau.Cor.Estimates.mat")).labpcexport;
cor_estimates_entropy_plaque = load(fullfile(results_dir,"Entropy.Plaque.Cor.Estimates.mat")).labpcexport;
cor_estimates_entropy_tau = load(fullfile(results_dir,"Entropy.Tau.Cor.Estimates.mat")).labpcexport;

% Set scale for the plots (must be the same scale)
overall_min = min([cor_estimates_TE_plaque(:); cor_estimates_TE_tau(:); cor_estimates_entropy_plaque(:); cor_estimates_entropy_tau(:)]); % find minimum to set as the lower bound
overall_max = max([cor_estimates_TE_plaque(:); cor_estimates_TE_tau(:); cor_estimates_entropy_plaque(:); cor_estimates_entropy_tau(:)]); % find the maximum to set as the upper bound

% Create brain plot for correlation between TE and plaque
close all;
data = cor_estimates_TE_plaque;
cmap = RdBu;

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[overall_min overall_max]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[overall_min overall_max]); % ensure clim matches clim on line 703
c = colorbar('WestOutside','fontsize',20); % color bar settings

% Save brain plot image
saveas(gcf,strcat(figures_dir,'Brain.Plot.TE.Plaque.Spearman.png'));
pause(2)



% Create brain plot for correlation between TE and tau
close all;
data = cor_estimates_TE_tau;
cmap = RdBu;

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[overall_min overall_max]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[overall_min overall_max]); % ensure clim matches clim on line 721
c = colorbar('WestOutside','fontsize',20); % color bar settings

% Save brain plot image
saveas(gcf,strcat(figures_dir,'Brain.Plot.TE.Tau.Spearman.png'));
pause(2)



% Create brain plot for correlation between entropy and plaque
close all;
data = cor_estimates_entropy_plaque;
cmap = RdBu;

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[overall_min overall_max]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[overall_min overall_max]); % ensure clim matches clim on line 739
c = colorbar('WestOutside','fontsize',20); % color bar settings

% Save brain plot image
saveas(gcf,strcat(figures_dir,'Brain.Plot.Entropy.Plaque.Spearman.png'));
pause(2)



% Create brain plot for correlation between entropy and tau
close all;
data = cor_estimates_entropy_tau;
cmap = RdBu;

img = display_atlas_blobs(data,'schaefer200x7','render',true,'backgroundimage',false,'colormap',cmap,'clim',[overall_min overall_max]); % can change clim if needed

figure('Position',[0 0 2000 2000]);
imshow(img);
set(gca,'colormap',colmap,'clim',[overall_min overall_max]); % ensure clim matches clim on line 757
c = colorbar('WestOutside','fontsize',20); % color bar settings

% Save brain plot image
saveas(gcf,strcat(figures_dir,'Brain.Plot.Entropy.Tau.Spearman.png'));
pause(2)



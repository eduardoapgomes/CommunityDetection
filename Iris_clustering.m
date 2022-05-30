%% Load Dataset to Workspace
clear;clc;rng('default');
load iris_dataset.mat
Dataset= irisInputs';
Target = irisTargets';
%% Preparing Data
% 1 - Normalization [0, 1]
Min    = min(Dataset);
Max    = max(Dataset);
Dataset = (Dataset - Min)./(Max-Min); 
%% Clustering
dm       = DAMICORE(Dataset,3);      % DAMICORE
KList    = 1:34;                     % Lista de Clusteres
Clusters = getclusterlist(dm,KList); % Extrair Resultados
%% Visualization
figure; 
subplot(5,1,1:4);dm.plotclusters();
subplot(5,1,5);dm.plotmodularity;grid on;
hold on;xlim([1,34]);
%% Cluster Evaluation
Eva1 = evalclusters(Dataset,Clusters,'CalinskiHarabasz');
Eva2 = evalclusters(Dataset,Clusters,'DaviesBouldin');
Eva3 = evalclusters(Dataset,Clusters,'silhouette');
%% Evaluation
figure;
subplot(311);Eva1.plot;grid on;hold on;xlim([1,34]);
subplot(312);Eva2.plot;grid on;hold on;xlim([1,34]);
subplot(313);Eva3.plot;grid on;hold on;xlim([1,34]);
%% Dataset Matrix
figure;
Clusters = getclusterlist(dm,3);   % (DAMICORE,List of Clusters) 
Leg = arrayfun(@(x)string(['C_{',int2str(x),'}']),Clusters);
[~,idx] = sort(Leg); 
gplotmatrix(Dataset(idx,:),...
    [],Leg(idx),linspecer(max(Clusters)),[],[],[],'grpbars');
grid on;
%% Table of Results
%Inputs
SepalLength = Dataset(:,1);
SepalWidth  = Dataset(:,2);
PetalLength = Dataset(:,3);
PetalWidth  = Dataset(:,4);   
%Outputs
Setosa      = Target(:,1);
Versicolor  = Target(:,2);
Virginica   = Target(:,3);
%Table
Inputs  = table(SepalLength,SepalWidth,PetalLength,PetalWidth,Clusters);
Outputs = table(Setosa,Versicolor,Virginica);
Table   = [Inputs,Outputs];
%% Parallel Plot:
figure; parallelplot(Table,'GroupVariable','Clusters');


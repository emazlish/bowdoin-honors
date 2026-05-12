% This is adapted code written by Sasha Kramer for use on HPLC pigments to
% be used with IFCB data
% EM 3/13/2026
%Map the directory where you will load your data:
cd /Users/emazlish/Library/CloudStorage/OneDrive-BowdoinCollege/F25/Honors/MatlabData/module5/

%Load your samples (formatted here as a .mat file):
load phytoSum_allsizes_noextrap_20260312.mat

tbl = phytoSum; % phytoSum for WSW samples, summer/fall Frac for >10 and >20 samples

%% create cluster table and clear pigment columns below threshold
IFCBcluster1 = tbl(:, {'sample','sumC_ugL','diatomC_ugL', 'dinoC_ugL', 'chrysoC_ugL', 'haptoC_ugL', 'cryptoC_ugL', 'prasinoC_ugL', 'eugC_ugL', 'cyanoC_ugL'});
% this does not include degradation pigment pheo a
labels = IFCBcluster1.Properties.VariableNames;
%Check percent of samples below detection (if >95% of values (EM version) are below detection, I
%remove the taxonomic group from the cluster analysis):
% for i = 2:10% columns with pigment info in them, i.e. chla through viola
%     belowThresh = IFCBcluster1(:,i) <= 0.001;
%     belowThresh = table2array(belowThresh);
%     pct(i) = 100*(sum(belowThresh)/size(IFCBcluster1, 1));
% end
% clear i
% idx2 = pct > 80; 

% put in new array
IFCBcluster2 = IFCBcluster1;
%IFCBcluster2(:,idx2) = [];
IFCBcluster2(:, 1:2) = []; % clear sample id and sumC_ugL
label2 = {'Diatom', 'Dinoflagellate', 'Chrysophyte', 'Haptophyte', 'Cryptophyte', 'Prasinophyte', 'Euglenoid'};
%label2(:, idx2) = [];
clear idx2
%% select only the rows you want for your cluster analysis by indexing
% original HPLC array
d1 = datetime('1/1/2023'); % start date
d2 = datetime('12/31/2025'); % end date
tf = (d1<tbl.HPLCdate & tbl.HPLCdate<d2); % select samples between dates
group = string(tbl.station) == 'CSC';
req = group ==1 & (tf == 1);
IFCBcluster2 = IFCBcluster2(req, :); % keep only those rows in rpigcluster2

%% Cluster group-based carbon, non normalized:
IFCBcluster2 = table2array(IFCBcluster2); % convert to array for clustering
D2 = pdist(IFCBcluster2','correlation'); %correlation is the method - you can vary that here
Z2 = linkage(D2,'ward'); %ward is the method - you can vary that here; I added 'euclidean' because I was getting a non-Euclidean distance warning

%Plot dendrogram of absolute pigment values:
figure;
h = dendrogram(Z2,'Labels',label2(1:end-1));
set(gca,'XTickLabelRotation',90,'fontsize',18)
set(h,'color','k','linewidth',2)
ylabel('Linkage Distance')
xlabel('Group')
title('Carbon by taxonomic group - non-normalized')
ax = gca;
ax.YGrid = 'on';
box on
clear ax h

%% Normalize to carbon and re-cluster:
normC = IFCBcluster2(:,1:end)./table2array(IFCBcluster1(req,2));
normC = table2array(normC);
normC(isnan(normC)) = 0;
normC(isinf(normC)) = 0;
normlabel = label2(1:end);

% match data with HPLC samples
IFCBtbl = phytoSum;
[idx, row] = ismember(Rpigcluster2.sample, IFCBtbl.sample(req));
normC = normC(row, :); % only keep 89 unique samples that also show up in pigment cluster for these dates
normC(:, 8) = []; % remove cyanos

% calculate correlation matrix
D3 = pdist(normC','correlation'); %correlation
Z3 = linkage(D3,'ward');
C3 = cophenet(Z3,D3);

%Plot dendrogram of normalized pigment values:
%figure();
h = dendrogram(Z3,'Labels',normlabel);
set(gca,'XTickLabelRotation',45,'fontsize',14)
set(h,'color','k','linewidth',2)
ylabel('Linkage distance')
xlabel('Group')
title('Carbon by taxonomic group - normalized to sumC')
ax = gca;
ax.YGrid = 'on';
box on

[c4, D] = cophenet(Z3, D3); % cophenetic correlation as done above with C3
[rho, pval] = corr(D3', D', "type", "Pearson"); % get rho (linear corr coeff) and pval from 'corr' function
hold on
hCircle = scatter(nan, nan, 40, 'filled', 'o', 'DisplayName', append('Pearson correlation coefficient = ', string(rho)));
hCircle = scatter(nan, nan, 40, 'filled', 'o', 'DisplayName', append('P-value = ', string(pval)));
numObs = length(normC);
legend({append('Pearson correlation coefficient = ', string(rho)), append('P-value = ', string(pval)), append('NumObservations = ', string(numObs))})
clear ax h D2 D3 deg other_data labels Z2 Z3 label2 label3
%% For further analysis: cluster samples - change "maxclust" based on your dendrogram results
%Need to also look at pigment concentrations in each cluster to check taxonomic meaning
C = clusterdata(normC,'distance','correlation','linkage','ward','maxclust',2);

% this will put the linkage distance cutoff on the dendrogram
%% p value based on Sasha Kramer, pers. comm. (1/2/25)
% use type 'Pearson' to calculate p value because 'cophenet' is doing a
% linear correlation, so we want the p value of the linear corr

[c4, D] = cophenet(Z3, D3); % cophenetic correlation as done above with C3
[rho, pval] = corr(D3', D', "type", "Pearson"); % get rho (linear corr coeff) and pval from 'corr' function


%% clustering by biovolume instead of carbon
IFCBcluster1 = tbl(:, {'sample','sumBiovolL','diatomVolL', 'dinoVolL', 'chrysoVolL', 'haptoVolL', 'cryptoVolL', 'prasinoVolL', 'eugVolL', 'cyanoVolL'});
% this does not include degradation pigment pheo a
labels = IFCBcluster1.Properties.VariableNames;
%Check percent of samples below detection (if >95% of values (EM version) are below detection, I
%remove the taxonomic group from the cluster analysis):
for i = 2:10% columns with pigment info in them, i.e. chla through viola
    belowThresh = IFCBcluster1(:,i) <= 0.001;
    belowThresh = table2array(belowThresh);
    percent(i) = 100*(sum(belowThresh)/size(IFCBcluster1, 1));
end
clear i
idx2 = percent > 95; 
% remove pigments below detection >80% of the time

% put in new array
IFCBcluster2 = IFCBcluster1;
IFCBcluster2(:,idx2) = [];
IFCBcluster2(:, 1:2) = []; % clear sample id and sumC_ugL

label2 = {'diatom', 'dino', 'chryso', 'hapto', 'crypto', 'prasino', 'eug', 'cyano'};
%label2(:, idx) = [];
clear idx2
%% select only the rows you want for your cluster analysis by indexing
% original HPLC array
d1 = datetime('1/1/2023'); % start date
d2 = datetime('12/31/2025'); % end date
tf = (d1<tbl.HPLCdate & tbl.HPLCdate<d2); % select samples between dates
group = string(tbl.station) == 'CSC';
req =  group == 1 &(tf == 1);
IFCBcluster2 = IFCBcluster2(req, :); % keep only those rows in rpigcluster2

%% Cluster group-based biovol:
IFCBcluster2 = table2array(IFCBcluster2); % convert to array for clustering
D2 = pdist(IFCBcluster2','correlation'); %correlation is the method - you can vary that here
Z2 = linkage(D2,'ward'); %ward is the method - you can vary that here; I added 'euclidean' because I was getting a non-Euclidean distance warning

%Plot dendrogram of absolute pigment values:
figure;
h = dendrogram(Z2,'Labels',label2([1:2, 4]));
set(gca,'XTickLabelRotation',90,'fontsize',18)
set(h,'color','k','linewidth',2)
ylabel('Linkage Distance')
xlabel('Group')
title('Biovolume by taxonomic group - non-normalized')
ax = gca;
ax.YGrid = 'on';
box on
clear ax h

%% Normalize to sumBiovol and re-cluster:
normBV = IFCBcluster2(:,1:end)./table2array(IFCBcluster1(req,2));
normBV(isnan(normBV)) = 0;
normBV(isinf(normBV)) = 0;
normlabel = label2%(1:end-1);

D3 = pdist(normBV','correlation'); %correlation
Z3 = linkage(D3,'ward');
C3 = cophenet(Z3,D3);

%Plot dendrogram of normalized pigment values:
figure();
h = dendrogram(Z3,'Labels',normlabel);
set(gca,'XTickLabelRotation',90,'fontsize',20)
set(h,'color','k','linewidth',2)
ylabel('Linkage distance')
xlabel('Group')
title('Biovolume by taxonomic group - normalized to sumBiovolL')
ax = gca;
ax.YGrid = 'on';
box on

[c4, D] = cophenet(Z3, D3); % cophenetic correlation as done above with C3
[rho, pval] = corr(D3', D', "type", "Pearson"); % get rho (linear corr coeff) and pval from 'corr' function
hold on
numObs = length(normBV);
legend({append('Pearson correlation coefficient = ', string(rho)), append('P-value = ', string(pval)), append('NumObservations = ', string(numObs))})

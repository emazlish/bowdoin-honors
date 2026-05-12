% This is adapted from code written by Sasha Kramer by Emma Mazlish
%Map the directory where you will load your data:
cd /Users/emazlish/Library/CloudStorage/OneDrive-BowdoinCollege/F25/Honors/MatlabData/module5/

%Load your samples (formatted here as a .mat file):
load phytoSum_allsizes_noextrap_20260312.mat HPLC

tbl = HPLC;
indC3 = (tbl.chlc3_ugl ~= 0 & tbl.chlc3_ugl <= 0.001); %Chlc3 0.001
tbl.chlc3_ugl(indC3) = 0; clear indC3

indC2 = (tbl.chlc2_ugl ~= 0 & tbl.chlc2_ugl <= 0.001); %Chlc2
tbl.chlc2_ugl(indC2)= 0; clear indC2

indC1 = (tbl.chlc1_ugl ~= 0 & tbl.chlc1_ugl <= 0.001); %Chlc1
tbl.chlc1_ugl(indC1) = 0; clear indC1

indVI = (tbl.viol_ugl ~= 0 & tbl.viol_ugl <= 0.001); %Viola
tbl.viol_ugl(indVI) = 0; clear indVI

indDD = (tbl.ddx_ugl ~= 0 & tbl.ddx_ugl <= 0.001); %Diadino
tbl.ddx_ugl(indDD) = 0; clear indDD

indDT = (tbl.dtx_ugl ~= 0 & tbl.dtx_ugl <= 0.001); %Diato
tbl.dtx_ugl(indDT) = 0; clear indDT

indAL = (tbl.allo_ugl ~= 0 & tbl.allo_ugl <= 0.001); %Allo
tbl.allo_ugl(indAL) = 0; clear indAL

indZE = (tbl.zea_ugl ~= 0 & tbl.zea_ugl <= 0.001); %Zea
tbl.zea_ugl(indZE) = 0; clear indZE

indLU = (tbl.lut_ugl ~= 0 & tbl.lut_ugl <= 0.001); %Lut
tbl.lut_ugl(indLU) = 0; clear indLU

%indPH = find(Global_RHPLC(:,24) ~= 0 & Global_RHPLC(:,24) <= 0.003); %Phide
%Global_RHPLC(indPH,24) = 0; clear indPH

indPE = (tbl.peri_ugl ~= 0 & tbl.peri_ugl <= 0.003); %Perid
tbl.peri_ugl(indPE) = 0; clear indPE

indBF = (tbl.but_ugl ~= 0 & tbl.but_ugl <= 0.002); %ButFuco
tbl.but_ugl(indBF) = 0; clear indBF

indFU = (tbl.fuco_ugl ~= 0 & tbl.fuco_ugl <= 0.002); %Fuco
tbl.fuco_ugl(indFU) = 0; clear indFU

indNE = find(tbl.neox_ugl ~= 0 & tbl.neox_ugl <= 0.002); %Neo
tbl.neox_ugl(indNE) = 0; clear indNE

indPR = find(tbl.pras_ugl ~= 0 & tbl.pras_ugl <= 0.002); %Pras
tbl.pras_ugl(indPR) = 0; clear indPR

indHF = (tbl.hex_ugl ~= 0 & tbl.hex_ugl <= 0.002); %HexFuco
tbl.hex_ugl(indHF) = 0; clear indHF

indCB = (tbl.chlb_ugl ~= 0 & tbl.chlb_ugl <= 0.003); %MVchlb
tbl.chlb_ugl(indCB) = 0; clear indCB

indDN = (tbl.dino_ugl ~= 0 & tbl.dino_ugl <= 0.001); %Dino
tbl.dino_ugl(indDN) = 0; clear indDN

indBE = (tbl.beta_ugl ~= 0 & tbl.beta_ugl <= 0.001); %Beta car
tbl.beta_ugl(indBE) = 0; clear indBE

indCHL = (tbl.chla_ugl ~= 0 & tbl.chla_ugl <= 0.001); %chla
tbl.chla_ugl(indCHL) = 0; clear indCHL

%% create cluster table and clear pigment columns below threshold
Rpigcluster1 = tbl(:, ["sample", "chla_ugl", "chlb_ugl", "beta_ugl", ...
    "but_ugl", "hex_ugl", "allo_ugl", "ddx_ugl", "dtx_ugl", "fuco_ugl", ...
    "peri_ugl", "zea_ugl", "chlc1_ugl", "chlc2_ugl", "chlc3_ugl", ...
    "lut_ugl", "neox_ugl", "dino_ugl", "pras_ugl", "viol_ugl"]);
% this does not include degradation pigment pheo a
labels = Rpigcluster1.Properties.VariableNames;
%Check percent of pigments below detection (if >80% samples are below detection, I
%remove the pigment from the cluster analysis):
for i = 2:20% columns with pigment info in them, i.e. chla through viola
    belowThresh = Rpigcluster1(:,i) <= 0.1;
    belowThresh = table2array(belowThresh);
    percent(i) = 100*(sum(belowThresh)/size(Rpigcluster1, 1));
end
clear i
idx = percent > 80; 
% remove pigments below detection >80% of the time

% put in new array
Rpigcluster2 = Rpigcluster1;
Rpigcluster2(:,idx) = [];
Rpigcluster2(:, 1) = []; % clear sample id
%Rpigcluster2(isnan(Rpigcluster2)) = 0;
%label2 = labels([2:17 19:20]); % remove sample label and dino
label2 = Rpigcluster2.Properties.VariableNames;
%label2 = {'chla', 'chlb', 'beta', 'but', 'hex', 'allo', 'ddx', 'dtx', 'fuco', 'peri', 'zea', 'chlc1', 'chlc2', 'chlc3', 'lut', 'neox', 'dino', 'pras', 'viol'};
%label2(:, idx) = [];
clear idx
%% select only the rows you want for your cluster analysis by indexing
% original HPLC array
d1 = datetime('1/1/2023'); % start date
d2 = datetime('12/31/2025'); % end date
tf = (d1<tbl.datetime & tbl.datetime<d2); % select samples between dates
group = string(tbl.station) == 'CSC' & (string(tbl.group) == 'A' | string(tbl.group) == 'B' | string(tbl.group) == 'CR' | string(tbl.group) == 'WSW');
req =  group == 1 & (tf == 1); %
Rpigcluster2 = Rpigcluster2(req, :); % keep only those rows in rpigcluster2
% for fall 2024 WSW, remove samples 34:35 and 39:40 because they are
% duplicates
%Rpigcluster2([12:13, 15:16], :) = [];
% then remove sample ID column
%Rpigcluster2(:, 1) = [];

matches = ismember(Rpigcluster2.sample, IFCBcluster2.sample);
%% Cluster pigments:
Rpigcluster2 = table2array(Rpigcluster2); % convert to array for clustering
D2 = pdist(Rpigcluster2','correlation'); %correlation is the method - you can vary that here
Z2 = linkage(D2,'ward'); %ward is the method - you can vary that here; I added 'euclidean' because I was getting a non-Euclidean distance warning

%Plot dendrogram of absolute pigment values:
figure(2),clf
h = dendrogram(Z2,'Labels',label2);
set(gca,'XTickLabelRotation',90,'fontsize',18)
set(h,'color','k','linewidth',2)
ylabel('Linkage Distance')
xlabel('Pigment')
title('Absolute pigment values')
ax = gca;
ax.YGrid = 'on';
box on
clear ax h

%% Normalize to chlorophyll-a and re-cluster:
normchl = Rpigcluster2(:,2:end)./Rpigcluster2(:,1);
normchl(isnan(normchl)) = 0;
normchl(isinf(normchl)) = 0;
normlabel = {'Chl b', 'Beta','But-fuco', 'Allo', 'Diadino', 'Fuco', 'Peri', 'Chl c1', 'Chl c2', 'Viol' };
%normlabel = label2(:, 2:end);

D3 = pdist(normchl','correlation'); %correlation
Z3 = linkage(D3,'ward');
C3 = cophenet(Z3,D3);

%Plot dendrogram of normalized pigment values:
figure();
h = dendrogram(Z3,'Labels',normlabel);
set(gca,'XTickLabelRotation',90,'fontsize',14)
set(h,'color','k','linewidth',2)
ylabel('Linkage distance')
xlabel('Pigment')
title('CSC WSW 2023-2025')
ax = gca;
ax.YGrid = 'on';
box on
hold on
[~, D] = cophenet(Z3, D3); % cophenetic correlation as done above with C3
[rho, pval] = corr(D3', D', "type", "Pearson"); % get rho (linear corr coeff) and pval from 'corr' function
numObs = length(normchl);
hold on
legend({append('Pearson correlation coefficient = ', string(rho)), append('P-value = ', string(pval)), append('NumObservations = ', string(numObs))})

clear ax h D2 D3 deg other_data labels Z2 Z3 label2 label3
%% For further analysis: cluster samples - change "maxclust" based on your dendrogram results
%Need to also look at pigment concentrations in each cluster to check taxonomic meaning
C = clusterdata(normchl,'distance','correlation','linkage','ward','maxclust',2);

% this will put the linkage distance cutoff on the dendrogram
%% p value based on Sasha Kramer's email (1/2/25)
% use type 'Pearson' to calculate p value because 'cophenet' is doing a
% linear correlation, so we want the p value of the linear corr

[c4, D] = cophenet(Z3, D3); % cophenetic correlation as done above with C3
[rho, pval] = corr(D3', D', "type", "Pearson"); % get rho (linear corr coeff) and pval from 'corr' function
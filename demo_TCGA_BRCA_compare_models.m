%% TCGA_BRCA_compare_models.m
% Read HiSeqV2 expression matrix + clinical; select top300 genes;
% split train/test; undersample training; run 3 models:
%  A) SVM RBF
%  B) PCA (k=50) + SVM RBF
%  C) Random Forest (TreeBagger)
% Compatible with MATLAB R2014a.

clc; clear; close all;
rng(1);

%% -------- USER SETTINGS --------
exprFile = 'HiSeqV2';   % tab-delimited expression matrix file (rows: genes, cols: samples; first col gene names)
clinFile = 'TCGA.BRCA.sampleMap_BRCA_clinicalMatrix'; % clinical file from Xena
topKGenes = 300;
pcaComponents = 50;     % for Model B; change if desired
rfNumTrees = 200;       % for Random Forest

outPrefix = 'TCGA_compare'; % prefix for saved PNGs/MATs

%% ===== 1. Read expression matrix =====
fprintf('Reading expression matrix: %s\n', exprFile);
fid = fopen(exprFile,'r');
if fid<0, error('Cannot open expression file %s', exprFile); end
headerLine = fgetl(fid);
headers = strsplit(headerLine, '\t');
numSamples = length(headers) - 1;

% read rest
fmt = ['%s' repmat('%f',1,numSamples)];
C = textscan(fid, fmt, 'Delimiter', '\t', 'TreatAsEmpty', {'NA',''} );
fclose(fid);

geneNames = C{1};
exprMatrix = cell2mat(C(2:end));  % genes x samples

fprintf('Expression size: %d genes x %d samples\n', size(exprMatrix,1), size(exprMatrix,2));

%% ===== 2. Read clinical (optional, we infer labels from sample codes) =====
fprintf('Reading clinical matrix (for metadata) : %s\n', clinFile);
fid = fopen(clinFile,'r');
if fid<0
    warning('Clinical file not found (%s). Continuing without clinical metadata.', clinFile);
    clinHeaders = {};
    clinicalSampleIDs = {};
else
    headerLine = fgetl(fid);
    clinHeaders = strsplit(headerLine, '\t');
    clinData = textscan(fid, repmat('%s',1,length(clinHeaders)), 'Delimiter', '\t', 'TreatAsEmpty', {'NA',''});
    fclose(fid);
    sampleCol = find(strcmpi(clinHeaders,'sampleID'),1);
    if ~isempty(sampleCol)
        clinicalSampleIDs = clinData{sampleCol};
    else
        clinicalSampleIDs = {};
    end
end

%% ===== 3. Build sample IDs and labels from expression header =====
exprSampleIDs = headers(2:end)'; % cell array of sample IDs (strings)
labels = zeros(length(exprSampleIDs),1);
for i=1:length(exprSampleIDs)
    s = exprSampleIDs{i};
    if length(s) >= 2
        code = s(end-1:end);
    else
        code = '00';
    end
    if strcmp(code,'01')
        labels(i) = 1;   % Tumor
    elseif strcmp(code,'11')
        labels(i) = 0;   % Normal
    else
        labels(i) = -1;  % Other
    end
end

validIdx = labels ~= -1;
exprMatrix = exprMatrix(:, validIdx);
labels = labels(validIdx);
exprSampleIDs = exprSampleIDs(validIdx);

fprintf('After filtering sample types: %d samples (Tumor=%d, Normal=%d)\n', ...
    length(labels), sum(labels==1), sum(labels==0));

%% ===== 4. Select topKGenes by variance =====
vars = var(exprMatrix,0,2);
[~, idxSort] = sort(vars,'descend');
topIdx = idxSort(1:min(topKGenes,length(idxSort)));
expr_reduced = exprMatrix(topIdx, :);    % genes x samples
selectedGeneNames = geneNames(topIdx);

% transpose to samples x features
X = expr_reduced';
Y = labels;

fprintf('Selected top %d genes ? feature matrix X size: %d samples x %d features\n', size(expr_reduced,1), size(X,1), size(X,2));

%% ===== 5. Train/test split =====
cv = cvpartition(Y, 'HoldOut', 0.2);
Xtrain = X(training(cv), :);
Ytrain = Y(training(cv));
Xtest  = X(test(cv), :);
Ytest  = Y(test(cv));

fprintf('Train samples: %d (Tumor=%d, Normal=%d)\n', length(Ytrain), sum(Ytrain==1), sum(Ytrain==0));
fprintf('Test samples:  %d (Tumor=%d, Normal=%d)\n', length(Ytest), sum(Ytest==1), sum(Ytest==0));

%% ===== 6. Undersample training set to balance classes =====
idxTum = find(Ytrain==1);
idxNorm = find(Ytrain==0);
nNorm = numel(idxNorm);
nTum  = numel(idxTum);
if nNorm==0 || nTum==0
    error('One of the classes has zero samples in training set. Check your labels.');
end
rng(1);
if nTum > nNorm
    idxTum_sub = randsample(idxTum, nNorm);
    balancedIdx = [idxNorm; idxTum_sub];
elseif nNorm > nTum
    idxNorm_sub = randsample(idxNorm, nTum);
    balancedIdx = [idxTum; idxNorm_sub];
else
    balancedIdx = [idxTum; idxNorm];
end

Xtrain_bal = Xtrain(balancedIdx, :);
Ytrain_bal = Ytrain(balancedIdx);

fprintf('Balanced training size: %d samples (Tumor=%d, Normal=%d)\n', length(Ytrain_bal), sum(Ytrain_bal==1), sum(Ytrain_bal==0));

%% small helper for confusion plotting (MATLAB2014a-safe)
% function plot_and_save_confusion(mat_counts, order, titleStr, fname)

%% ===== Model A: SVM RBF (no PCA) =====
fprintf('\n--- Model A: SVM (RBF) ---\n');
MdlA = fitcsvm(Xtrain_bal, Ytrain_bal, 'KernelFunction', 'rbf', 'Standardize', true, 'KernelScale','auto');
YpredA = predict(MdlA, Xtest);
accA = mean(YpredA == Ytest);
C_A = confusionmat(Ytest, YpredA);
% compute precision/recall/F1 for class 1 (Tumor)
orderA = unique([0;1]);
TP = C_A(orderA==1, orderA==1);
FP = C_A(orderA==0, orderA==1);
FN = C_A(orderA==1, orderA==0);
precisionA = TP / (TP + FP + eps);
recallA    = TP / (TP + FN + eps);
F1_A       = 2*precisionA*recallA/(precisionA+recallA+eps);

fprintf('Accuracy (A) = %.2f%%\n', accA*100);
fprintf('Precision = %.3f, Recall = %.3f, F1 = %.3f\n', precisionA, recallA, F1_A);
% save confusion images (counts + normalized)
plot_and_save_confusion(C_A, orderA, 'SVM RBF - counts', [outPrefix '_SVM_counts.png']);
C_A_norm = bsxfun(@rdivide, C_A, sum(C_A,2)+eps);  % row-normalized
plot_and_save_confusion(C_A_norm, orderA, 'SVM RBF - row-normalized', [outPrefix '_SVM_norm.png']);

%% ===== Model B: PCA + SVM RBF =====
fprintf('\n--- Model B: PCA (%d comps) + SVM (RBF) ---\n', pcaComponents);
% PCA on training data (centered)
[coeff, score_train, ~, ~, ~] = pca(Xtrain);
k = min(pcaComponents, size(score_train,2));
Xtrain_pca = score_train(:,1:k);
% project Xtest
Xtest_centered = Xtest - repmat(mean(Xtrain), size(Xtest,1),1);
Xtest_pca = Xtest_centered * coeff(:,1:k);

% balance on PCA space (use the same balanced indices from original balancedIdx)
% Need to recompute balanced indices mapping relative to training rows
train_inds = find(training(cv)); % indices in original X that went to train
train_inds = train_inds(:);
balanced_train_inds = train_inds(balancedIdx); % absolute indices in X
% Now create Xtrain_pca_bal, Ytrain_pca_bal
% find positions of balanced_train_inds within training set ordering => they are balancedIdx already
Xtrain_pca_bal = Xtrain_pca(balancedIdx, :);
Ytrain_pca_bal = Ytrain_bal;

MdlB = fitcsvm(Xtrain_pca_bal, Ytrain_pca_bal, 'KernelFunction','rbf','Standardize',true,'KernelScale','auto');
YpredB = predict(MdlB, Xtest_pca);
accB = mean(YpredB == Ytest);
C_B = confusionmat(Ytest, YpredB);
TP = C_B(orderA==1, orderA==1);
FP = C_B(orderA==0, orderA==1);
FN = C_B(orderA==1, orderA==0);
precisionB = TP / (TP + FP + eps);
recallB    = TP / (TP + FN + eps);
F1_B       = 2*precisionB*recallB/(precisionB+recallB+eps);

fprintf('Accuracy (B) = %.2f%%\n', accB*100);
fprintf('Precision = %.3f, Recall = %.3f, F1 = %.3f\n', precisionB, recallB, F1_B);
plot_and_save_confusion(C_B, orderA, 'PCA+SVM - counts', [outPrefix '_PCA_SVM_counts.png']);
C_B_norm = bsxfun(@rdivide, C_B, sum(C_B,2)+eps);
plot_and_save_confusion(C_B_norm, orderA, 'PCA+SVM - row-normalized', [outPrefix '_PCA_SVM_norm.png']);

%% ===== Model C (fixed) - TreeBagger (no OOBPrediction) =====
fprintf('\n--- Model C: Random Forest (TreeBagger) - variant A ---\n');

Xtrain_rf = Xtrain_bal;
Ytrain_rf = Ytrain_bal; % numeric 0/1

% convert labels to cell-array of strings (TreeBagger handles classification with string labels reliably)
Ytrain_rf_cell = cellstr(num2str(Ytrain_rf));

% Train
TB = TreeBagger(rfNumTrees, Xtrain_rf, Ytrain_rf_cell, 'Method','classification');

% Predict on Xtest (predict returns cell array of label-strings)
[YpredC_cell, scores] = predict(TB, Xtest);
YpredC = str2double(YpredC_cell); % convert back to numeric 0/1

% Evaluate
accC = mean(YpredC == Ytest);
C_C = confusionmat(Ytest, YpredC);
TP = C_C(orderA==1, orderA==1);
FP = C_C(orderA==0, orderA==1);
FN = C_C(orderA==1, orderA==0);
precisionC = TP / (TP + FP + eps);
recallC    = TP / (TP + FN + eps);
F1_C       = 2*precisionC*recallC/(precisionC+recallC+eps);

fprintf('Accuracy (C - TreeBagger) = %.2f%%\n', accC*100);
fprintf('Precision = %.3f, Recall = %.3f, F1 = %.3f\n', precisionC, recallC, F1_C);

plot_and_save_confusion(C_C, orderA, 'RandomForest - counts', [outPrefix '_RF_counts.png']);
C_C_norm = bsxfun(@rdivide, C_C, sum(C_C,2)+eps);
plot_and_save_confusion(C_C_norm, orderA, 'RandomForest - row-normalized', [outPrefix '_RF_norm.png']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ===== Model C: Random Forest (TreeBagger) =====
% fprintf('\n--- Model C: Random Forest (TreeBagger) ---\n');
% % TreeBagger expects predictors as samples x features, and responses categorical
% Ytrain_rf = Ytrain_bal; % balanced labels
% Xtrain_rf = Xtrain_bal;
% 
% % Use TreeBagger
% TB = TreeBagger(rfNumTrees, Xtrain_rf, Ytrain_rf, 'OOBPrediction','On', 'Method','classification');
% [YpredC_cell, scores] = predict(TB, Xtest); % returns cell array of strings '0'/'1'
% YpredC = str2double(YpredC_cell);
% 
% accC = mean(YpredC == Ytest);
% C_C = confusionmat(Ytest, YpredC);
% TP = C_C(orderA==1, orderA==1);
% FP = C_C(orderA==0, orderA==1);
% FN = C_C(orderA==1, orderA==0);
% precisionC = TP / (TP + FP + eps);
% recallC    = TP / (TP + FN + eps);
% F1_C       = 2*precisionC*recallC/(precisionC+recallC+eps);
% 
% fprintf('Accuracy (C) = %.2f%%\n', accC*100);
% fprintf('Precision = %.3f, Recall = %.3f, F1 = %.3f\n', precisionC, recallC, F1_C);
% plot_and_save_confusion(C_C, orderA, 'RandomForest - counts', [outPrefix '_RF_counts.png']);
% C_C_norm = bsxfun(@rdivide, C_C, sum(C_C,2)+eps);
% plot_and_save_confusion(C_C_norm, orderA, 'RandomForest - row-normalized', [outPrefix '_RF_norm.png']);
% 

%% ===== Save results =====
results.SVM = struct('acc',accA,'precision',precisionA,'recall',recallA,'F1',F1_A,'confusion',C_A);
results.PCA_SVM = struct('acc',accB,'precision',precisionB,'recall',recallB,'F1',F1_B,'confusion',C_B);
results.RF = struct('acc',accC,'precision',precisionC,'recall',recallC,'F1',F1_C,'confusion',C_C);
save([outPrefix '_results.mat'], 'results', 'selectedGeneNames');

fprintf('\nAll models done. Results saved to %s_results.mat and confusion PNGs.\n', outPrefix);

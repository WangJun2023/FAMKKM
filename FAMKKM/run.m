close all; clear all; clc
warning off;

ResSavePath = 'FAMKKM/result/';
dataPath = 'Datasets/Kernel-Datasets/';
datasetName = {'Pollen', 'caltech101_nTrain10_48_Kmatrix','caltech101_nTrain20_48_Kmatrix','caltech101_nTrain30_48_Kmatrix','flower102_Kmatrix',...
    'ALOI_100_Kmatrix'};

for dataIndex = 3 : 6
    dataName = [dataPath datasetName{dataIndex} '.mat'];
    load(dataName);
    
    KH = kcenter(KH);
    KH = knorm(KH);
    numClust = length(unique(Y));
    
    ResBest = zeros(1, 8);
    ResStd = zeros(1, 8);
    
    % parameters setting
    r1 = [0.01, 0.1, 1];
    r2 = [0.01, 0.1, 1];

    acc = zeros(length(r1), length(r2));
    nmi = zeros(length(r1), length(r2));
    purity = zeros(length(r1), length(r2));
    idx = 1;
    for r1Index = 1 : length(r1)
        r1Temp = r1(r1Index);
        for r2Index = 1 : length(r2)
            r2Temp = r2(r2Index);
            fprintf('Please wait a few minutes\n');
            disp(['Dataset: ', datasetName{dataIndex}, ...
                ', --r1--: ', num2str(r1Temp), ', --r2--: ', num2str(r2Temp)]);
            tic;
            [res,H] = demo_main(KH, numClust, Y, r1Temp, r2Temp);
            Runtime(idx) = toc;
            disp(['runtime: ', num2str(Runtime(idx))]);
            idx = idx + 1;
            tempResBest(1, : ) = res(1, : );
            tempResStd(1, : ) = res(2, : );
            acc(r1Index, r2Index) = tempResBest(1, 7);
            nmi(r1Index, r2Index) = tempResBest(1, 4);
            purity(r1Index, r2Index) = tempResBest(1, 8);
            for tempIndex = 1 : 8
                if tempResBest(1, tempIndex) > ResBest(1, tempIndex)
                    ResBest(1, tempIndex) = tempResBest(1, tempIndex);
                    ResStd(1, tempIndex) = tempResStd(1, tempIndex);
                end
            end
        end
    end
    aRuntime = mean(Runtime);
    resFile2 = [ResSavePath datasetName{dataIndex}, '-ACC=', num2str(ResBest(1, 7)), '.mat'];
    save(resFile2, 'ResBest', 'ResStd', 'acc', 'nmi', 'purity', 'aRuntime');
end
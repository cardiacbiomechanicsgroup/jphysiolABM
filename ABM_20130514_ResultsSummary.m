clear all;

% subdirectory to search for results files
FileList = {'ABM_20130514PAR_Biax_Circle_AD_R3_Ws0.16667_42days.mat' 'ABM_20130514PAR_UniaxC_Circle_AD_R3_Ws0.16667_42days.mat'};

% name of file to save the following analysis
svnm = 'ABM_20130314_ResultsSummary.mat';

% number of files in results directory
NumFiles = 2;

lst = 1:1:NumFiles;                       % specify order of analyzing and plotting results files
leg = {'Biax, Circle' 'Uniax, Circle'};   % specify simulation parameters for each simulation
fmt = {'-r' ':b'};                        % specify plotting line style for each simulation

T = (0:1:42)';  % simulation time (days)

% arrays for storing time course results
CellTimeMVA = zeros([43 NumFiles]);  % mean angle of cells within infarct region versus time
CellTimeMVL = zeros([43 NumFiles]);  % mean vector length of cells within infarct region versus time
CellTimeFRC = zeros([43 NumFiles]);  % area fraction of cells within infarct region versus time
CollTimeMVA = zeros([43 NumFiles]);  % mean angle of collagen fibers within infarct region versus time
CollTimeMVL = zeros([43 NumFiles]);  % mean vector length of collagen fibers within infarct region versus time
CollTimeFRC = zeros([43 NumFiles]);  % area fraction of collagen fibers within infarct region versus time
FbrnTimeFRC = zeros([43 NumFiles]);  % area fraction of non-collagen fibers ("fibrin") within infarct region versus time

% arrays for storing cell and collagen fiber orientation histograms 3 weeks
% after infarction
CellHist3WK = zeros([36 NumFiles]);
CollHist3WK = zeros([36 NumFiles]);



for ifile = 1:1:NumFiles

    % load results file of first simulation
    load(FileList{ifile}, '-mat');


    
    NumTimeSamples = floor((NumTimeSteps+1)/(SamplingInterval/TimeStep)) + 1;



    % calculate fiber patch positions
    [x, y] = meshgrid(xmin:FiberPatchSpacing:xmax, ymin:FiberPatchSpacing:ymax);
    NumColFiberPatches = numel(x);
    [R, C] = size(x);

    ColFiberPosX = reshape(x,NumColFiberPatches,1);
    ColFiberPosY = reshape(y,NumColFiberPatches,1);

    ColFiberZone = zeros([NumColFiberPatches 1]);

    ColFiberDistances = zeros([NumColFiberPatches 1]);

    AveSpan = 3;
    CFPX = zeros([floor(R/AveSpan) floor(C/AveSpan)]);
    CFPY = zeros([floor(R/AveSpan) floor(C/AveSpan)]);
    for i = 1:AveSpan:R-rem(R,AveSpan)
        for j = 1:AveSpan:C-rem(C,AveSpan)
            CFPX((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan) = mean2(x(i:i+AveSpan-1,j:j+AveSpan-1));
            CFPY((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan) = mean2(y(i:i+AveSpan-1,j:j+AveSpan-1));
        end
    end

    clear x y;


    
    % results file contains NumSims = 4 repeats of each simulation.  arrays
    % collect the data for each repeat.
    CellMVA = zeros([NumTimeSamples NumSims]);
    CellMVL = zeros([NumTimeSamples NumSims]);
    CellFRC = zeros([NumTimeSamples NumSims]);
    CollMVA = zeros([NumTimeSamples NumSims]);
    CollMVL = zeros([NumTimeSamples NumSims]);
    CollFRC = zeros([NumTimeSamples NumSims]);
    FbrnFRC = zeros([NumTimeSamples NumSims]);
    Cell3wkHist = zeros([NumFiberDirBins NumSims]);
    Coll3wkHist = zeros([NumFiberDirBins NumSims]);

    for i = 1:1:NumSims

        CellMVA(:,i) = zMeanWoundCellAngle{i}(:,1);
        CellMVL(:,i) = zMeanWoundCellLength{i}(:,1);
        CellFRC(:,i) = zNCWound{i}(:,1);
        CollMVA(:,i) = zMeanWoundColFiberAngle{i}(:,1);
        CollMVL(:,i) = zMeanWoundColFiberLength{i}(:,1);
        CollFRC(:,i) = zNCFWound{i}(:,1);
        FbrnFRC(:,i) = zNFFWound{i}(:,1);
        Cell3wkHist(:,i) = mean(zWoundCellAngleDistribution{i}(:,20:24,1)./repmat(sum(zWoundCellAngleDistribution{i}(:,20:24,1)), NumFiberDirBins, 1), 2);
        Coll3wkHist(:,i) = mean(zWoundColFiberAngleDistribution{i}(:,20:24,1)./repmat(sum(zWoundColFiberAngleDistribution{i}(:,20:24,1)), NumFiberDirBins, 1), 2);

    end

    
    
    % calculate mean vector components
    CellMVX = CellMVL.*cos(2*CellMVA*pi/180);
    CellMVY = CellMVL.*sin(2*CellMVA*pi/180);
    CollMVX = CollMVL.*cos(2*CollMVA*pi/180);
    CollMVY = CollMVL.*sin(2*CollMVA*pi/180);

    
    
    % average the data for all of the repeat simulations.  note, mean
    % angles and mean vector lengths cannot be averaged directly.  the mean
    % vector components must be averaged and then converted back to a mean
    % angle and mean vector length.  standard deviations of the mean angle
    % and mean vector length are derived from the standard error
    % propagation formula.
    % means
    CellMVXm = mean(CellMVX,2);
    CellMVYm = mean(CellMVY,2);
    CellFRCm = mean(CellFRC,2)/WoundArea*CellArea;
    CollMVXm = mean(CollMVX,2);
    CollMVYm = mean(CollMVY,2);
    CollFRCm = mean(CollFRC,2)/WoundArea/FiberDensity;
    FbrnFRCm = mean(FbrnFRC,2)/WoundArea/FiberDensity;
    Cell3wkHistm = mean(Cell3wkHist,2);
    Coll3wkHistm = mean(Coll3wkHist,2);

    % standard deviations
    CellMVXs = std(CellMVX,0,2);
    CellMVYs = std(CellMVY,0,2);
    CellFRCs = std(CellFRC,0,2)/WoundArea*CellArea;
    CollMVXs = std(CollMVX,0,2);
    CollMVYs = std(CollMVY,0,2);
    CollFRCs = std(CollFRC,0,2)/WoundArea/FiberDensity;
    FbrnFRCs = std(FbrnFRC,0,2)/WoundArea/FiberDensity;
    Cell3wkHists = std(Cell3wkHist,0,2);
    Coll3wkHists = std(Coll3wkHist,0,2);

    % means
    CellMVAm = 180/pi*1/2*atan2(CellMVYm,CellMVXm);
    CellMVLm = sqrt(CellMVYm.^2 + CellMVXm.^2);
    CollMVAm = 180/pi*1/2*atan2(CollMVYm,CollMVXm);
    CollMVLm = sqrt(CollMVYm.^2 + CollMVXm.^2);

    % standard deviations
    CellMVAs = sqrt(CellMVYs.^2.*(180/pi*1/2*CellMVXm./CellMVLm.^2).^2 + CellMVXs.^2.*(180/pi*1/2*-1*CellMVYm./CellMVLm.^2).^2);
    CellMVLs = sqrt(CellMVYs.^2.*(CellMVYm./CellMVLm).^2 + CellMVXs.^2.*(CellMVXm./CellMVLm).^2);
    CollMVAs = sqrt(CollMVYs.^2.*(180/pi*1/2*CollMVXm./CollMVLm.^2).^2 + CollMVXs.^2.*(180/pi*1/2*-1*CollMVYm./CollMVLm.^2).^2);
    CollMVLs = sqrt(CollMVYs.^2.*(CollMVYm./CollMVLm).^2 + CollMVXs.^2.*(CollMVXm./CollMVLm).^2);

    
    
    % store averaged data
    CellTimeMVA(1:NumTimeSamples,ifile) = CellMVAm;
    CellTimeMVL(1:NumTimeSamples,ifile) = CellMVLm;
    CellTimeFRC(1:NumTimeSamples,ifile) = CellFRCm;
    CollTimeMVA(1:NumTimeSamples,ifile) = CollMVAm;
    CollTimeMVL(1:NumTimeSamples,ifile) = CollMVLm;
    CollTimeFRC(1:NumTimeSamples,ifile) = CollFRCm;
    FbrnTimeFRC(1:NumTimeSamples,ifile) = FbrnFRCm;
    
    CellHist3WK(:,ifile) = Cell3wkHistm;
    CollHist3WK(:,ifile) = Coll3wkHistm;
    
    
    
    % plot the collagen field for the first repeat simulation
    figure(1);
    if ifile == 1
    clf;
    end
    if NumFiles > 3
    subplot(2,ceil(NumFiles/2),ifile);
    else
    subplot(1,NumFiles,ifile);
    end
    hold on;
    contourf(CFPX, CFPY, zCFND{1}(:,:,22)/FiberPatchSpacing^2/FiberDensity, 10, 'LineStyle', 'none');
    colormap((3*jet + white)/4);
    caxis([0 0.3]);
    plot(WoundPerimX, WoundPerimY, '-r', 'LineWidth', 2);
    quiver(CFPX, CFPY, zCFMVX{1}(:,:,22), zCFMVY{1}(:,:,22), 0.67, 'k', 'ShowArrowHead', 'off', 'LineWidth', 1);
    hold off;
    title([leg{ifile} '   Wk3 Collagen'], 'FontName', 'Arial', 'FontSize', 14);
    ylabel('Y (\mum)', 'FontName', 'Arial', 'FontSize', 14);
    xlabel('X (\mum)', 'FontName', 'Arial', 'FontSize', 14);
    axis([xmin xmax ymin ymax]);
    set(gca, 'DataAspectRatio', [1 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);
    
    
    
    % plot the cell field for the first repeat simulation
    figure(2);
    if ifile == 1
    clf;
    end
    if NumFiles > 3
    subplot(2,ceil(NumFiles/2),ifile);
    else
    subplot(1,NumFiles,ifile);
    end
    hold on;
    plot(WoundPerimX, WoundPerimY, '-b', 'LineWidth', 2);
    plot(zCPX{1}{22}, zCPY{1}{22}, 'or', 'MarkerSize', 10, 'LineWidth', 1);
    quiver(zCPX{1}{22}, zCPY{1}{22}, zCDX{1}{22}, zCDY{1}{22}, 0.25, 'r', 'ShowArrowHead', 'off', 'LineWidth', 1);
    hold off;
    title([leg{ifile} '   Wk3 Cells'], 'FontName', 'Arial', 'FontSize', 14);
    ylabel('Y (\mum)', 'FontName', 'Arial', 'FontSize', 14);
    xlabel('X (\mum)', 'FontName', 'Arial', 'FontSize', 14);
    axis([xmin xmax ymin ymax]);
    set(gca, 'DataAspectRatio', [1 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);
            
end



% plot the time course data
figure(3);
clf;

subplot(2,3,1);
hold on;
for i = 1:1:numel(lst)
plot(T, CollTimeMVL(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
% title(' ', 'FontName', 'Arial', 'FontSize', 14);
ylabel('MVL_C_F', 'FontName', 'Arial', 'FontSize', 14);
xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 14);
axis([0 50 0 1]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);

subplot(2,3,2);
hold on;
for i = 1:1:numel(lst)
plot(T, CollTimeMVA(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
% title(' ', 'FontName', 'Arial', 'FontSize', 14);
ylabel('MVA_C_F', 'FontName', 'Arial', 'FontSize', 14);
xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 14);
axis([0 50 -90 90]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);

subplot(2,3,3);
hold on;
for i = 1:1:numel(lst)
plot(T, CollTimeFRC(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
% title(' ', 'FontName', 'Arial', 'FontSize', 14);
ylabel('F_C_F', 'FontName', 'Arial', 'FontSize', 14);
xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 14);
legend(leg, 'Location', 'NorthWest');
axis([0 50 0 1]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);

subplot(2,3,4);
hold on;
for i = 1:1:numel(lst)
plot(T, CellTimeMVL(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
% title(' ', 'FontName', 'Arial', 'FontSize', 14);
ylabel('MVL_C_E_L_L', 'FontName', 'Arial', 'FontSize', 14);
xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 14);
axis([0 50 0 1]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);

subplot(2,3,5);
hold on;
for i = 1:1:numel(lst)
plot(T, CellTimeMVA(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
% title(' ', 'FontName', 'Arial', 'FontSize', 14);
ylabel('MVA_C_E_L_L', 'FontName', 'Arial', 'FontSize', 14);
xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 14);
axis([0 50 -90 90]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);

subplot(2,3,6);
hold on;
for i = 1:1:numel(lst)
plot(T, CellTimeFRC(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
% title(' ', 'FontName', 'Arial', 'FontSize', 14);
ylabel('F_C_E_L_L', 'FontName', 'Arial', 'FontSize', 14);
xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 14);
axis([0 50 0 1]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);



% plot the collagen and cell orientation histograms
figure(4);
clf;
subplot(1,2,1);
hold on;
for i = 1:1:numel(lst)
plot(FiberDirBins, CollHist3WK(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
title('3wk Collagen', 'FontName', 'Arial', 'FontSize', 14);
ylabel('Frequency', 'FontName', 'Arial', 'FontSize', 14);
xlabel('\theta_C_O_L_L (deg)', 'FontName', 'Arial', 'FontSize', 14);
axis([-90 90 0 0.1]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);

subplot(1,2,2);
hold on;
for i = 1:1:numel(lst)
plot(FiberDirBins, CellHist3WK(:,lst(i)), fmt{i}, 'LineWidth', 2);
end
hold off;
title('3wk Cells', 'FontName', 'Arial', 'FontSize', 14);
ylabel('Frequency', 'FontName', 'Arial', 'FontSize', 14);
xlabel('\theta_C_E_L_L (deg)', 'FontName', 'Arial', 'FontSize', 14);
legend(leg, 'Location', 'NorthWest');
axis([-90 90 0 0.1]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);






save(svnm, '-mat');



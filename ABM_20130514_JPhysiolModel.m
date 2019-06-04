clear all;

% Model setup

% Mechanical environments that can be implemented. 'UniaxC' implements a
% strain field in the infarct region where Exx (circumferential) is
% +0.05 and Eyy (longitudinal) is zero.  'UniaxL' implements a strain
% field in the infarct region where Eyy is +0.05 and Exx is zero. 'Biax'
% implements a strain field in the infarct region where Exx=Eyy=0.05.  In
% all three cases the strains in the remote tissue are Exx=Eyy=-0.05.
mechs = {'UniaxC' 'UniaxL' 'Biax'};

% Fiber deposition processes that can be implemented. 'AD' implements
% aligned deposition, where cells deposit collagen fibers that are aligned
% with the current cell orientation. 'RD' implments random deposition.
dep = {'AD' 'RD'};

% Infarct shapes that can be implemented. Elliptical shapes can have an
% aspect ratio of 2 or 5 and be oriented circumferentially (Horizontal) or
% longitudinally (Vertical).  'Ellipse2_Vertical' and 'Ellipse2_Horizontal'
% model the experimental elliptical cryoinfarcts.
shapes = {'Circle' 'Ellipse2_Vertical' 'Ellipse2_Horizontal' 'Ellipse5_Vertical' 'Ellipse5_Horizontal'};


% Number of times to repeat the simulation
NumSims = 4;


% Length and migration rate scaling factors
% 113um is model circular cryoinfarct radius.  2828um is experimental
% circular cryoinfarct radius.  500um is experimental infarct
% half-thickness.
AlphaGeometric = 113/2828;
AlphaInfiltration = 113/500;


% Persistence times for calculation of persistence tuning factor (later)
Pmin = 0.25;     % lower bound on persistence time (hr)
Pnocue = 1.2;    % persistence time in the absence of external directional cues (hr)


% Time discretization.  'SamplingInterval' specifies when certain model
% outputs should be saved.
TimeStep = 2*Pmin;                             % hr
NumDays = 42;                                  % d
NumTimeSteps = round(NumDays*24/TimeStep);     % number of iterations
Times = (0:TimeStep:NumTimeSteps*TimeStep)';   % hr
SamplingInterval = 24;                         % hr
SamplingInterval2 = 24*7;                      % hr


% Theta space discretization for various orientation-based calculations.
Directions = (-179:1:180)';
DirectionLine = (-180:1:180)';
DirectionBins = (-180+15/2:15:180-15/2)';
DirBins360 = (-177.5:5:177.5)';
FiberDirBins = (-90+5/2:5:90-5/2)';
NumFiberDirBins = numel(FiberDirBins);


% Fibroblast parameters
CellRadius = 5;     % um
CellArea = pi*CellRadius^2;    % um2

CellSpacing = 20;   % um

CellDensityMax = 0.0085;    % cells/um2 (maximum concentration of cells in the model with overlap prohibited)
CellAreaFracMax = CellDensityMax*CellArea;    % cell area / total area

CellSpeedMin = AlphaInfiltration*1;     % um/hr
CellSpeedMax = AlphaInfiltration*10;    % um/hr

CellCycleDurationMax = 240;   % hr 
CellCycleDurationMin = 24;    % hr
ApoptosisThreshold = 240;     % hr

% X and Y coordinates of cell border and unit outward normal vector
% components at these coordinates. 'Local' means coordinates are relative
% to the cell centroid.
CellPerimLocalX = CellRadius*cos(DirectionLine*pi/180);
CellPerimLocalY = CellRadius*sin(DirectionLine*pi/180);
CellNormalX = cos(DirectionLine*pi/180);
CellNormalY = sin(DirectionLine*pi/180);

       
% Fiber parameters (Rotation, Generation, and Degradation rate
% coefficients)
FiberPatchSpacing = CellRadius/2;

FiberDensity = 178.25;    % fibers/um2/7um thickness

kColFiberRotMax = 3;  % 0.75;           % deg/hr
kColFiberRotMin = 0.1*kColFiberRotMax;           % deg/hr

% (FiberDensity*CellArea) converts generation rate from [fiber area fraction/hr]
% to [fibers/cell/hr].  (1.18/CellAreaFracMax) adjusts generation rate to
% account for the fact that cells do not uniformly cover the infarct
% area and do not deposit collagen uniformly over the infarct area
% (which is assumed in the simple ODE for collagen accumulation).
kColFiberGenMax = 0.0007269*FiberDensity*CellArea*1.18/CellAreaFracMax;            % fibers/cell/hr
kColFiberGenMin = 0.000007269*FiberDensity*CellArea*1.18/CellAreaFracMax;          % fibers/cell/hr

% (1/CellAreaFracMax) adjusts the degradation rate to account for the fact
% that cells do not uniformly cover the infarct area.
kColFiberDegMin = 0.00025152/CellAreaFracMax;           % unitless scaling factor
kColFiberDegMax = 0.0025152/CellAreaFracMax;             % unitless scaling factor

kFibFiberDegMax = 0.0025152/CellAreaFracMax;     % hr-1
kFibFiberDegMin = 0.00025152/CellAreaFracMax;    % hr-1

  
% Chemokine parameters
ConcMax = 1;  % need to know these from the concentration field that is imported later
ConcMin = 0;

% Maximum possible concentration gradient magnitude, for the prescribed
% chemokine diffusion, generation, and degradation parameters, for
% normalizing the chemical guidance cue vector.
kcdeg = 0.001;         % 1/s
kcgen = 0.01;          % nM/s
Dc = 100;              % um2/s
L = 2828;              % um
lam = sqrt(kcdeg/Dc);  % 1/um

DimlessMaxConcGradMag = lam*L*sinh(lam*L)/(cosh(lam*L) + sinh(lam*L));

MaxConcGradMag = DimlessMaxConcGradMag/(AlphaGeometric*L);    % 1/um    note, concentration is still normalized


% Mechanics parameters
% Strains in remote tissue.
EpsilonTissueXX = -0.05;
EpsilonTissueYY = -0.05;
EpsilonTissueXY = 0;


% Functional border zone as fraction of average infarct diameter.
% Prescribes the width of the transition of mechanics-related parameters
% (strains and degradation rate coefficient) from the infarct values to the
% remote values across the infarct border.
FBZFracXX = 0.08;
FBZFracYY = 0.08;
FBZFracXY = 0.08;


% Maximum possible strain anisotropy, based on cryoinfarction data, for
% normalizing the mechanical guidance cue vector.
eii = 0.05;
ejj = 0;
eij = 0;


% Guidance cue weight factors
Wc = 1/3;  % cell persistence
Wcf = 1/6;  % collagen fibers
ws = 1/6*[1 3.5 7.5];  % controlling strain cue weight factor in a for loop
Wcg = 1/6;  % chemokine gradient


%Guidance cue normalization factors
Mc = CellSpeedMax;
Mcf = 1;
Ms = sqrt(eii^2 - 2*eii*ejj + ejj^2 + 4*eij^2)/4;
Mcg = CellRadius*MaxConcGradMag/2;


% Persistence tuning factor.  Enforces persistent random walk behavior.
RhoTuningFactor = 2*Wc*Pmin/(Pnocue - Pmin);




for imechs = [3 1]  % 1:1:numel(mechs)               % 1. UniaxC,  2. UniaxL,  3. Biax

    for idep = 1  % 1:1:numel(dep)               % 1. AD,  2. RD

        % Specify aligned deposition or random deposition
        deposition = dep{idep};

        for ishapes = 1 % 1:1:numel(shapes)      % 1. Circle,  2. Ellipse2_Vertical,  3. Ellipse2_Horizontal,  4. Ellipse5_Vertical,  5. Ellipse5_Horizontal

            for iws = 1 % 1:1:numel(ws)          % 1/6*[1 3.5 7.5]

                % Specify strain cue weight factor
                Ws = ws(iws);



                savename = ['ABM_20130514_' mechs{imechs} '_' shapes{ishapes} '_' deposition '_R' num2str(kColFiberRotMax) '_Ws' num2str(Ws) '_' num2str(NumDays) 'days'];


                % Collect output from each repeat of the simulation
                % in these cell arrays
                zT = cell([1 NumSims]);  % time of each saved time step
                zNC = cell([1 NumSims]);  % total number of cells at each saved time step
                zNCF = cell([1 NumSims]);  % total number of collagen fibers at each saved time step
                zCWI = cell([1 NumSims]);  % list of indices of cells located within the infarct wound region at each saved time step
                zCPX = cell([1 NumSims]);  % cell x-coordinates at each saved time step
                zCPY = cell([1 NumSims]);  % cell y-coordinates at each saved time step
                zCDX = cell([1 NumSims]);  % cell x-orientations at each saved time step
                zCDY = cell([1 NumSims]);  % cell y-orientations at each saved time step
                zCFWI = cell([1 NumSims]);  % list of indices of collagen fiber patches (elements) located within the infarct wound region at each saved time step
                zCFMVX = cell([1 NumSims]);  % locally averaged collagen fiber mean vector x-component at each saved time step
                zCFMVY = cell([1 NumSims]);  % locally averaged collagen fiber mean vector y-component at each saved time step
                zCFND = cell([1 NumSims]);  % locally averaged mean number of collagen fibers per patch (element) at each saved time step
                zFFND = cell([1 NumSims]);  % locally averaged mean number of fibrin fibers per patch (element) at each saved time step
                zNCWound = cell([1 NumSims]);  % number of cells within the infarct wound region at each saved time step
                zMeanWoundCellAngle = cell([1 NumSims]);  % mean orientation of cells within the infarct wound region at each saved time step
                zMeanWoundCellLength = cell([1 NumSims]);  % mean vector length of cell orientation vectors within the infarct wound region at each saved time step
                zWoundCellAngleDistribution = cell([1 NumSims]);  % cell orientation histogram of cells within the infarct wound region at each saved time step
                zNCFWound = cell([1 NumSims]);  % number of collagen fibers within the infarct wound region at each saved time step
                zMeanWoundColFiberAngle = cell([1 NumSims]);  % mean orientation of collagen fibers within the infarct wound region at each saved time step
                zMeanWoundColFiberLength = cell([1 NumSims]);  % mean vector length of collagen orientation vectors within the infarct wound region at each saved time step
                zWoundColFiberAngleDistribution = cell([1 NumSims]);  % collagen fiber orientation histogram of collagen fibers within the infarct wound region at each saved time step
                zNFFWound = cell([1 NumSims]);  % number of fibrin fibers within the infarct wound region at each saved time step
                zT2 = cell([1 NumSims]);  % time of each saved time step
                zColPX = cell([1 NumSims]);  % every collagen fiber patch (element) x-coordinates
                zColPY = cell([1 NumSims]);  % every collagen fiber patch (element) y-coordinates
                zColDX = cell([1 NumSims]);  % every collagen fiber patch (element) mean vector x-component
                zColDY = cell([1 NumSims]);  % every collagen fiber patch (element) mean vector y-component


                % Compute wound strain tensor
                if strcmp(mechs{imechs}, 'UniaxL')

                    EpsilonWoundXX = 0;
                    EpsilonWoundYY = 0.05;
                    EpsilonWoundXY = 0;

                elseif strcmp(mechs{imechs}, 'UniaxC')

                    EpsilonWoundXX = 0.05;
                    EpsilonWoundYY = 0;
                    EpsilonWoundXY = 0;

                else

                    EpsilonWoundXX = 0.05;
                    EpsilonWoundYY = 0.05;
                    EpsilonWoundXY = 0;

                end




                % Load chemokine concentration profile and wound dimensions, 
                % and define subregions within and around wound
                load(['ABM_ConcField_' shapes{ishapes} '.mat'], '-mat', 'X', 'Y', 'Conc', 'WoundRadiusX', 'WoundRadiusY');

                WoundPerimX = (WoundRadiusX*WoundRadiusY./sqrt((WoundRadiusX*sin(DirectionLine*pi/180)).^2 + (WoundRadiusY*cos(DirectionLine*pi/180)).^2)).*cos(DirectionLine*pi/180);
                WoundPerimY = (WoundRadiusX*WoundRadiusY./sqrt((WoundRadiusX*sin(DirectionLine*pi/180)).^2 + (WoundRadiusY*cos(DirectionLine*pi/180)).^2)).*sin(DirectionLine*pi/180);

                WoundArea = pi*WoundRadiusX*WoundRadiusY;

                % functional border zone width (um)
                StdvFBZXX = FBZFracXX*sqrt(4*WoundRadiusX*WoundRadiusY);
                StdvFBZYY = FBZFracYY*sqrt(4*WoundRadiusX*WoundRadiusY);
                StdvFBZXY = FBZFracXY*sqrt(4*WoundRadiusX*WoundRadiusY);

                % dimensions of 2D space (um)
                buffer = 100;

                xmax = ceil((WoundRadiusX + buffer)/10)*10;
                xmin = -xmax;
                ymax = ceil((WoundRadiusY + buffer)/10)*10;
                ymin = -ymax;
                limit = ceil((max([WoundRadiusX WoundRadiusY]) + buffer)/10)*10;
                axislimits = round(1.1*[-limit limit -limit limit]);

                clear buffer limit;





%                 matlabpool open;
%                 parfor repeat = 1:NumSims
                for repeat = 1:1:NumSims


                    % initialize cells
                    NumCells = round(((xmax - xmin)/CellSpacing)*((ymax - ymin)/CellSpacing));

                    CellPosX = rand([NumCells 1])*(xmax - xmin - 2*CellRadius) + xmin + CellRadius;
                    CellPosY = rand([NumCells 1])*(ymax - ymin - 2*CellRadius) + ymin + CellRadius;

                    CellZone = zeros([NumCells 1]);

                    Dummy = 2*pi*rand([NumCells 1]);
                    CellDirX = cos(Dummy);
                    CellDirY = sin(Dummy);
    %                 clear Dummy;

                    % eliminate cells from wound region
                    CellPosR = sqrt(CellPosX.^2 + CellPosY.^2);
                    CellPosT = atan2(CellPosY, CellPosX);
                    WoundPerimR = WoundRadiusX*WoundRadiusY./sqrt((WoundRadiusX*sin(CellPosT)).^2 + (WoundRadiusY*cos(CellPosT)).^2);

                    CellPosX = CellPosX(CellPosR >= WoundPerimR);
                    CellPosY = CellPosY(CellPosR >= WoundPerimR);
                    CellDirX = CellDirX(CellPosR >= WoundPerimR);
                    CellDirY = CellDirY(CellPosR >= WoundPerimR);

                    WPR = WoundPerimR(CellPosR >= WoundPerimR);
                    CPR = CellPosR(CellPosR >= WoundPerimR);

                    NumCells = numel(CellPosX);

    %                 clear CellPosR CellPosT WoundPerimR;


                    % initialize fibroblast dynamic parameters
                    Concentration = interp2(X, Y, Conc, CellPosX, CellPosY);

                    CellCycleDuration = (Concentration - ConcMin)/(ConcMax - ConcMin)*(CellCycleDurationMin - CellCycleDurationMax) + CellCycleDurationMax;

                    CellSpeed = (Concentration - ConcMin)/(ConcMax - ConcMin)*(CellSpeedMax - CellSpeedMin) + CellSpeedMin;

                    CellCycleClock = round(rand([NumCells 1]).*CellCycleDuration./TimeStep).*TimeStep;

                    ApoptosisClock = CellCycleClock;

                    kColFiberRot = (Concentration - ConcMin)/(ConcMax - ConcMin)*(kColFiberRotMax - kColFiberRotMin) + kColFiberRotMin;

                    kColFiberGen = (Concentration - ConcMin)/(ConcMax - ConcMin)*(kColFiberGenMax - kColFiberGenMin) + kColFiberGenMin;

                    kColFiberDeg = (Concentration - ConcMin)/(ConcMax - ConcMin)*(kColFiberDegMax - kColFiberDegMin) + kColFiberDegMin;
                    
                    kFibFiberDeg = (Concentration - ConcMin)/(ConcMax - ConcMin)*(kFibFiberDegMax - kFibFiberDegMin) + kFibFiberDegMin;


    %                 clear Concentration CPR WPR;


                    % initialize fibers
                    [x, y] = meshgrid(xmin:FiberPatchSpacing:xmax, ymin:FiberPatchSpacing:ymax);
                    NumColFiberPatches = numel(x);
                    [R, C] = size(x);

                    ColFiberPosX = reshape(x,NumColFiberPatches,1);
                    ColFiberPosY = reshape(y,NumColFiberPatches,1);

                    ColFiberZone = zeros([NumColFiberPatches 1]);

                    ColFiberDistances = zeros([NumColFiberPatches 1]);

                    % locally averaged collagen fiber patch (element) coordinates
                    AveSpan = 3;
                    CFPX = zeros([floor(R/AveSpan) floor(C/AveSpan)]);
                    CFPY = zeros([floor(R/AveSpan) floor(C/AveSpan)]);
                    for i = 1:AveSpan:R-rem(R,AveSpan)
                        for j = 1:AveSpan:C-rem(C,AveSpan)
                            CFPX((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan) = mean2(x(i:i+AveSpan-1,j:j+AveSpan-1));
                            CFPY((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan) = mean2(y(i:i+AveSpan-1,j:j+AveSpan-1));
                        end
                    end

    %                 clear x y;

                    % initialize fiber orientation histograms for
                    % each patch (element)

                    % for isotropic fiber distribution
                    %{
                    ColFiberDistribution = ones([NumColFiberPatches NumFiberDirBins]);
                    FibrinFiberDistribution = 25*ones([NumColFiberPatches NumFiberDirBins]);
                    %}

                    % for anisotropic fiber distribution
                    %{
                    ColFiberDistribution = zeros([NumColFiberPatches NumFiberDirBins]);
                    ColFiberDistribution(:,18:19) = 18;
                    FibrinFiberDistribution = zeros([NumColFiberPatches NumFiberDirBins]);
                    FibrinFiberDistribution(:,18:19) = 450;
                    %}

                    % for native fiber distribution

                    NativeFiberDistribution = [0.000043 0.000079 0.000174 0.000377 0.000782 0.001542 0.002885 0.005124 0.008638 0.013821 0.020988 0.03025 0.041382 0.053729 0.06621 0.07744 0.085965 0.090574 0.090574 0.085965 0.07744 0.06621 0.053729 0.041382 0.03025 0.020988 0.013821 0.008638 0.005124 0.002885 0.001542 0.000782 0.000377 0.000174 0.000079 0.000043];
                    ColFiberDistribution = 36*repmat(NativeFiberDistribution, NumColFiberPatches, 1);
                    FibrinFiberDistribution = 25*36*repmat(NativeFiberDistribution, NumColFiberPatches, 1);
    %                 clear NativeFiberDistribution;


                    % find indices of collagen fiber patches
                    % (elements) located within the infarct wound
                    % region
                    ColFiberPosR = sqrt(ColFiberPosX.^2 + ColFiberPosY.^2);
                    ColFiberPosT = atan2(ColFiberPosY, ColFiberPosX);
                    WoundPerimR = WoundRadiusX*WoundRadiusY./sqrt((WoundRadiusX*sin(ColFiberPosT)).^2 + (WoundRadiusY*cos(ColFiberPosT)).^2);

                    WoundColFiberIndices = find(ColFiberPosR <= WoundPerimR);

                    CFWI = WoundColFiberIndices;

    %                 clear ColFiberPosR ColFiberPosT WoundPerimR InnerPerimR OuterPerimR;


                    % for isotropic fiber distribution in wound
                    %{
                    ColFiberDistribution(WoundColFiberIndices,:) = 1;
                    FibrinFiberDistribution(WoundColFiberIndices,:) = 25;
                    %}


                    % collagen orientation statistics for each fiber patch
                    % mean angle, mean vector length, and mean
                    % vector components
                    DirectionBinsTiled = repmat(FiberDirBins', NumColFiberPatches, 1);
                    Cbar = sum(ColFiberDistribution.*cos(2*DirectionBinsTiled*pi/180), 2)./sum(ColFiberDistribution, 2);
                    Sbar = sum(ColFiberDistribution.*sin(2*DirectionBinsTiled*pi/180), 2)./sum(ColFiberDistribution, 2);

                    ColFiberMA = 180/pi*1/2*atan2(Sbar,Cbar);
                    ColFiberMVL = sqrt(Cbar.^2 + Sbar.^2);
                    ColFiberMVX = ColFiberMVL.*cos(ColFiberMA*pi/180);
                    ColFiberMVY = ColFiberMVL.*sin(ColFiberMA*pi/180);

    %                 clear DirectionBinsTiled Cbar Sbar;


                    % zone assignment for faster processing of
                    % "neighborhood" queries.  zones are larger
                    % than fiber patches (elements).  assign fiber
                    % patches and cells to the zones in which they
                    % reside.  generate connectivity list
                    % ("AdjacentZones") for each zone.
                    [x, y] = meshgrid(xmin:10:xmax-10, ymin:10:ymax-10);
                    MinX = reshape(x,numel(x),1);
                    MinY = reshape(y,numel(y),1);

                    [x, y] = meshgrid(xmin+10:10:xmax, ymin+10:10:ymax);
                    MaxX = reshape(x,numel(x),1);
                    MaxY = reshape(y,numel(y),1);

                    [r, c] = size(x);

    %                 clear x y;

                    AdjacentZones = cell([numel(MinX) 1]);

                    ZoneColFiberIndices = cell([numel(MinX) 1]);

                    for i = 1:1:numel(MinX)

                        ColFiberZone(ColFiberPosX >= MinX(i) & ColFiberPosX < MaxX(i) & ColFiberPosY >= MinY(i) & ColFiberPosY < MaxY(i)) = i;
                        ZoneColFiberIndices{i} = find(ColFiberPosX >= MinX(i) & ColFiberPosX < MaxX(i) & ColFiberPosY >= MinY(i) & ColFiberPosY < MaxY(i));
                        CellZone(CellPosX >= MinX(i) & CellPosX < MaxX(i) & CellPosY >= MinY(i) & CellPosY < MaxY(i)) = i;

                        if i == 1
                            AdjacentZones{i} = [i+1 i+r i+r+1];
                        elseif i > 1 && i < r
                            AdjacentZones{i} = [i-1 i+1 i+r-1 i+r i+r+1];
                        elseif i == r
                            AdjacentZones{i} = [i-1 i+r-1 i+r];
                        elseif mod(i-1,r) == 0 && i < numel(MinX)-r
                            AdjacentZones{i} = [i-r i-r+1 i+1 i+r i+r+1];
                        elseif mod(i,r) == 0 && i < numel(MinX)
                            AdjacentZones{i} = [i-r-1 i-r i-1 i+r-1 i+r];
                        elseif i == numel(MinX)-r+1
                            AdjacentZones{i} = [i-r i-r+1 i+1];
                        elseif i == numel(MinX)
                            AdjacentZones{i} = [i-r-1 i-r i-1];
                        elseif i > numel(MinX)-r+1 && i < numel(MinX)
                            AdjacentZones{i} = [i-r-1 i-r i-r+1 i-1 i+1];
                        else
                            AdjacentZones{i} = [i-r-1 i-r i-r+1 i-1 i+1 i+r-1 i+r i+r+1];
                        end

                    end




                    % initialize variables for data storage:
                    % migration and alignment metrics.  all of
                    % these variables correspond to those with the
                    % "z" prefix defined above, but these are just
                    % storing the information for this single
                    % repeat of the simulation.
                    NumTimeSamples = floor((NumTimeSteps+1)/(SamplingInterval/TimeStep)) + 1;

                    T = zeros([1 NumTimeSamples]);
                    NC = zeros([NumTimeSamples 1]);
                    NCF = zeros([NumTimeSamples 1]);
                    CPX = cell([1 NumTimeSamples]);
                    CPY = cell([1 NumTimeSamples]);
                    CDX = cell([1 NumTimeSamples]);
                    CDY = cell([1 NumTimeSamples]);
                    CWI = cell([NumTimeSamples 1]);
                    CFMVX = zeros([size(CFPX) NumTimeSamples]);
                    CFMVY = zeros([size(CFPX) NumTimeSamples]);
                    CFND = zeros([size(CFPX) NumTimeSamples]);
                    FFND = zeros([size(CFPX) NumTimeSamples]);
                    NCWound = zeros([NumTimeSamples 1]);
                    MeanWoundCellAngle = zeros([NumTimeSamples 1]);
                    MeanWoundCellLength = zeros([NumTimeSamples 1]);
                    WoundCellAngleDistribution = zeros([numel(FiberDirBins) NumTimeSamples 1]);
                    NCFWound = zeros([NumTimeSamples 1]);
                    MeanWoundColFiberAngle = zeros([NumTimeSamples 1]);
                    MeanWoundColFiberLength = zeros([NumTimeSamples 1]);
                    WoundColFiberAngleDistribution = zeros([numel(FiberDirBins) NumTimeSamples 1]);
                    NFFWound = zeros([NumTimeSamples 1]);

                    NumTimeSamples2 = floor((NumTimeSteps+1)/(SamplingInterval2/TimeStep)) + 1;

                    T2 = zeros([1 NumTimeSamples2]);
                    ColPX = zeros([NumColFiberPatches NumTimeSamples2]);
                    ColPY = zeros([NumColFiberPatches NumTimeSamples2]);
                    ColDX = zeros([NumColFiberPatches NumTimeSamples2]);
                    ColDY = zeros([NumColFiberPatches NumTimeSamples2]);




                    % store initial state

                    ntp = 1;

                    T(ntp) = 0;

                    NC(ntp) = NumCells;
                    NCF(ntp) = sum(sum(ColFiberDistribution));

                    CPX{ntp} = CellPosX;
                    CPY{ntp} = CellPosY;
                    CDX{ntp} = CellDirX;
                    CDY{ntp} = CellDirY;


                    % locally averaged collagen fiber orientation
                    % statistics and number density
                    CFMAtemp = reshape(ColFiberMA, R, C);
                    CFMVLtemp = reshape(ColFiberMVL, R, C);
                    CFMVXtemp = reshape(ColFiberMVX, R, C);
                    CFMVYtemp = reshape(ColFiberMVY, R, C);
                    CFNDtemp = reshape(sum(ColFiberDistribution,2), R, C);

                    FFNDtemp = reshape(sum(FibrinFiberDistribution,2), R, C);

                    for i = 1:AveSpan:R-rem(R,AveSpan)
                        for j = 1:AveSpan:C-rem(C,AveSpan)

                            Ctemp = sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1).*CFMVLtemp(i:i+AveSpan-1,j:j+AveSpan-1).*cos(2*CFMAtemp(i:i+AveSpan-1,j:j+AveSpan-1)*pi/180)))/sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1)));
                            Stemp = sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1).*CFMVLtemp(i:i+AveSpan-1,j:j+AveSpan-1).*sin(2*CFMAtemp(i:i+AveSpan-1,j:j+AveSpan-1)*pi/180)))/sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1)));
                            Atemp = 1/2*atan2(Stemp, Ctemp);
                            Mtemp = sqrt(Stemp^2 + Ctemp^2);

                            CFMVX((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = Mtemp*cos(Atemp);
                            CFMVY((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = Mtemp*sin(Atemp);
                            CFND((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = mean2(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1));

                            FFND((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = mean2(FFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1));

    %                         clear Ctemp Stemp Atemp Mtemp;

                        end
                    end

    %                 clear CFMVXtemp CFMVYtemp CFNDtemp FFNDtemp CFMAtemp CFMVLtemp;


                    % find indices of cells located within the
                    % various infarct wound subregions
                    CellPosR = sqrt(CellPosX.^2 + CellPosY.^2);
                    CellPosT = atan2(CellPosY, CellPosX);
                    WoundPerimR = WoundRadiusX*WoundRadiusY./sqrt((WoundRadiusX*sin(CellPosT)).^2 + (WoundRadiusY*cos(CellPosT)).^2);

                    WoundCellIndices = find(CellPosR <= WoundPerimR);

                    CWI{ntp,1} = WoundCellIndices;


                    % infarct wound average cell orientation
                    % measures and number of cells.  repeat for all
                    % wound subregions.
                    if ~isempty(WoundCellIndices)

                        NCWound(ntp,1) = numel(WoundCellIndices);

                        WoundCellDirX = CellDirX(WoundCellIndices);
                        WoundCellDirY = CellDirY(WoundCellIndices);
                        WoundCellAngle = atan2(WoundCellDirY, WoundCellDirX);

                        MeanWoundCellDirX = mean(cos(2*WoundCellAngle));
                        MeanWoundCellDirY = mean(sin(2*WoundCellAngle));

                        MeanWoundCellAngle(ntp,1) = 180/pi*1/2*atan2(MeanWoundCellDirY, MeanWoundCellDirX);
                        MeanWoundCellLength(ntp,1) = sqrt(MeanWoundCellDirX^2 + MeanWoundCellDirY^2);

                        WoundCellAngle(WoundCellAngle > pi/2 & WoundCellAngle <= pi) = WoundCellAngle(WoundCellAngle > pi/2 & WoundCellAngle <= pi) - pi;
                        WoundCellAngle(WoundCellAngle > -pi & WoundCellAngle <= -pi/2) = WoundCellAngle(WoundCellAngle > -pi & WoundCellAngle <= -pi/2) + pi;
                        WoundCellAngleDistribution(:,ntp,1) = hist(WoundCellAngle*180/pi, FiberDirBins);

                    else

                        NCWound(ntp,1) = 0;
                        MeanWoundCellAngle(ntp,1) = 0;
                        MeanWoundCellLength(ntp,1) = 0;
                        WoundCellAngleDistribution(:,ntp,1) = 0;

                    end

    %                 clear CellPosR CellPosT WoundPerimR WoundCellDirX WoundCellDirY WoundCellAngle MeanWoundCellDirX MeanWoundCellDirY;


                    % infarct wound average collagen fiber
                    % orientation measures and number of fibers.
                    % repeat for all wound subregions.
                    NCFWound(ntp,1) = sum(sum(ColFiberDistribution(WoundColFiberIndices,:)));
                    NFFWound(ntp,1) = sum(sum(FibrinFiberDistribution(WoundColFiberIndices,:)));

                    if NCFWound(ntp,1) > 0

                        WoundColFiberAngleDistribution(:,ntp,1) = (sum(ColFiberDistribution(WoundColFiberIndices,:), 1))';

                        MeanWoundColFiberDirX = sum(WoundColFiberAngleDistribution(:,ntp,1).*cos(2.*FiberDirBins.*pi./180))./sum(WoundColFiberAngleDistribution(:,ntp,1));
                        MeanWoundColFiberDirY = sum(WoundColFiberAngleDistribution(:,ntp,1).*sin(2.*FiberDirBins.*pi./180))./sum(WoundColFiberAngleDistribution(:,ntp,1));

                        MeanWoundColFiberAngle(ntp,1) = 180/pi*1/2*atan2(MeanWoundColFiberDirY, MeanWoundColFiberDirX);
                        MeanWoundColFiberLength(ntp,1) = sqrt(MeanWoundColFiberDirX^2 + MeanWoundColFiberDirY^2);

                    else

                        MeanWoundColFiberAngle(ntp,1) = 0;
                        MeanWoundColFiberLength(ntp,1) = 0;
                        WoundColFiberAngleDistribution(:,ntp,1) = 0;

                    end

    %                 clear MeanWoundColFiberDirX MeanWoundColFiberDirY;

                    ntp2 = 1;

                    T2(ntp2) = 0;

                    % store every collagen fiber patch position and
                    % orientation mean vector components
                    ColPX(:,ntp2) = ColFiberPosX;
                    ColPY(:,ntp2) = ColFiberPosY;
                    ColDX(:,ntp2) = ColFiberMVX;
                    ColDY(:,ntp2) = ColFiberMVY;

                    % plot
                    %{
                    figure(1)
                    clf;
                    hold on;
                    quiver(ColFiberPosX, ColFiberPosY, ColFiberMVX, ColFiberMVY, 1, 'k', 'LineWidth', 1);
                    quiver(CellPosX, CellPosY, CellDirX, CellDirY, 0.3, 'r', 'LineWidth', 3);
                    plot(CellPosX, CellPosY, 'or', 'MarkerSize', 10, 'LineWidth', 2);
                    plot(WoundPerimX, WoundPerimY, '-b', 'LineWidth', 4);
                    axis(axislimits);
                    axis square;
                    grid on;
                    title('Fibers and Cells', 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                    ylabel('Y Position (um)', 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                    xlabel('X Position (um)', 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                    set(gca, 'LineWidth', 2, 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                    hold off;
                    pause;
                    %}







                    % start simulation
                    for t = 2:1:NumTimeSteps+1

    %                     tic;

                        % clear indices of cells flagged for apoptosis
                        ApoptosisIndex = [];

                        % randomize order of looping through cells
                        CurrNumCells = NumCells;
                        RandLine = rand([CurrNumCells 1]);
                        [~, RandCellIndex] = sort(RandLine);
    %                     clear RandLine;


                        for ic = 1:1:CurrNumCells

                            i = RandCellIndex(ic);





                            % cell death
                            if ApoptosisClock(i) >= ApoptosisThreshold

                                ApoptosisIndex = [ApoptosisIndex i];

                            end





                            % cell division
                            if CellCycleClock(i) >= CellCycleDuration(i)

                                % calculate available directions
                                AvailableDirections = ones(size(Directions));

                                CellDistances = sqrt((CellPosX - CellPosX(i)).^2 + (CellPosY - CellPosY(i)).^2);

                                NearbyCells = find(CellDistances > 0 & CellDistances < 4*CellRadius);

                                if ~isempty(NearbyCells)

                                    for j = 1:1:numel(NearbyCells)

                                        Xnew = CellPosX(i) + 2*CellRadius*cos(Directions*pi/180);
                                        Ynew = CellPosY(i) + 2*CellRadius*sin(Directions*pi/180);

                                        DistanceNewToNearby = sqrt((CellPosX(NearbyCells(j)) - Xnew).^2 + (CellPosY(NearbyCells(j)) - Ynew).^2);

                                        AvailableDirections(DistanceNewToNearby < 2*CellRadius) = 0;

                                    end

                                end

                                PossibleDirections = Directions(AvailableDirections == 1);

    %                             clear CellDistances NearbyCells Xnew Ynew DistanceNewToNearby;                            


                                if sum(AvailableDirections) > 0

                                    % create a new cell in a randomly selected available location
                                    NumCells = NumCells + 1;

                                    RandomlySelectedDir = PossibleDirections(randi([1,numel(PossibleDirections)],1));

                                    CellPosX(NumCells) = CellPosX(i) + 2*CellRadius*cos(RandomlySelectedDir*pi/180);
                                    CellPosY(NumCells) = CellPosY(i) + 2*CellRadius*sin(RandomlySelectedDir*pi/180);

    %                                 clear RandomlySelectedDir;

                                    CellDirX(NumCells) = CellDirX(i);
                                    CellDirY(NumCells) = CellDirY(i);
                                    CellSpeed(NumCells) = CellSpeed(i);
                                    CellCycleClock(NumCells) = round(0.1*CellCycleDurationMin*randn/TimeStep)*TimeStep;
                                    CellCycleDuration(NumCells) = CellCycleDuration(i);
                                    ApoptosisClock(NumCells) = CellCycleClock(NumCells);
                                    CellZone(NumCells) = CellZone(i);

                                    CellCycleClock(i) = round(0.1*CellCycleDurationMin*randn/TimeStep)*TimeStep;
                                    ApoptosisClock(i) = CellCycleClock(i);

                                else

                                    % no space available for division, cell goes quiescent for
                                    % a time step.
                                    CellCycleClock(i) = round(0.1*CellCycleDurationMin*randn/TimeStep)*TimeStep;
                                    ApoptosisClock(i) = ApoptosisClock(i) + TimeStep;

                                end

    %                             clear AvailableDirections PossibleDirections;






                            % cell reorientation and migration, and collagen
                            % reorientation, degradation, and deposition
                            else

                                % increment clocks
                                CellCycleClock(i) = CellCycleClock(i) + TimeStep;
                                ApoptosisClock(i) = ApoptosisClock(i) + TimeStep;


                                % compute global coordinates around cell
                                CellPerimGlobalX = CellPerimLocalX + CellPosX(i);
                                CellPerimGlobalY = CellPerimLocalY + CellPosY(i);


                                % compute direction and magnitude of cell persistence
                                CX = CellDirX(i);
                                CY = CellDirY(i);
                                CM = CellSpeed(i);


                                % compute cell cycle duration, cell speed, and
                                % collagen turnover parameters
                                Concentration = interp2(X, Y, Conc, CellPerimGlobalX, CellPerimGlobalY);

                                ConcMean = mean(Concentration);

                                CellCycleDuration(i) = (ConcMean - ConcMin)/(ConcMax - ConcMin)*(CellCycleDurationMin - CellCycleDurationMax) + CellCycleDurationMax;

                                CellSpeed(i) = (ConcMean - ConcMin)/(ConcMax - ConcMin)*(CellSpeedMax - CellSpeedMin) + CellSpeedMin;

                                kColFiberRot(i) = (ConcMean - ConcMin)/(ConcMax - ConcMin)*(kColFiberRotMax - kColFiberRotMin) + kColFiberRotMin;

                                kColFiberGen(i) = (ConcMean - ConcMin)/(ConcMax - ConcMin)*(kColFiberGenMax - kColFiberGenMin) + kColFiberGenMin;

                                kColFiberDeg(i) = (ConcMean - ConcMin)/(ConcMax - ConcMin)*(kColFiberDegMax - kColFiberDegMin) + kColFiberDegMin;
                                
                                kFibFiberDeg(i) = (ConcMean - ConcMin)/(ConcMax - ConcMin)*(kFibFiberDegMax - kFibFiberDegMin) + kFibFiberDegMin;


                                % compute chemical guidance cue vector
                                MeanGradientX = 1/2/pi*trapz(DirectionLine*pi/180, Concentration.*cos(DirectionLine*pi/180));
                                MeanGradientY = 1/2/pi*trapz(DirectionLine*pi/180, Concentration.*sin(DirectionLine*pi/180));

                                CGM = sqrt(MeanGradientX^2 + MeanGradientY^2);  % chemokine gradient cue vector magnitude
                                CGX = MeanGradientX/CGM;  % chemokine gradient cue unit vector x-component
                                CGY = MeanGradientY/CGM;  % chemokine gradient cue unit vector y-component

    %                             clear ConcMean CellPosR CellPosT WoundPerimR MeanGradientX MeanGradientY;


                                % compute strain guidance cue vector
                                CellPerimGlobalR = sqrt(CellPerimGlobalX.^2 + CellPerimGlobalY.^2);
                                CellPerimGlobalT = atan2(CellPerimGlobalY, CellPerimGlobalX);
                                WoundBorderR = WoundRadiusX*WoundRadiusY./sqrt((WoundRadiusX*sin(CellPerimGlobalT)).^2 + (WoundRadiusY*cos(CellPerimGlobalT)).^2);

                                StrainXX = EpsilonWoundXX + (EpsilonTissueXX - EpsilonWoundXX)*1/2*erfc(-(CellPerimGlobalR-WoundBorderR)/StdvFBZXX/sqrt(2));
                                StrainYY = EpsilonWoundYY + (EpsilonTissueYY - EpsilonWoundYY)*1/2*erfc(-(CellPerimGlobalR-WoundBorderR)/StdvFBZYY/sqrt(2));
                                StrainXY = EpsilonWoundXY + (EpsilonTissueXY - EpsilonWoundXY)*1/2*erfc(-(CellPerimGlobalR-WoundBorderR)/StdvFBZXY/sqrt(2));

                                NormalStrain = (StrainXX.*CellNormalX + StrainXY.*CellNormalY).*CellNormalX + (StrainXY.*CellNormalX + StrainYY.*CellNormalY).*CellNormalY;

                                MeanStrainX = 1/2/pi*trapz(DirectionLine*pi/180, NormalStrain.*cos(2*DirectionLine*pi/180));
                                MeanStrainY = 1/2/pi*trapz(DirectionLine*pi/180, NormalStrain.*sin(2*DirectionLine*pi/180));
                                MeanStrainAngle = 1/2*atan2(MeanStrainY, MeanStrainX);

                                SX = cos(MeanStrainAngle);  % strain cue unit vector x-component
                                SY = sin(MeanStrainAngle);  % strain cue unit vector y-component
                                SM = sqrt(MeanStrainX^2 + MeanStrainY^2);  % strain cue vector magnitude

    %                             clear Concentration StrainXX StrainYY StrainXY NormalStrain MeanStrainX MeanStrainY MeanStrainAngle CellPerimGlobalX CellPerimGlobalY CellPerimGlobalR CellPerimGlobalT WoundBorderR;


                                % compute collagen fiber guidance cue vector
                                % (this actually includes collagen and fibrin fibers)
                                ColFiberDistances(:) = 1000*CellRadius;

                                AdjacentColFiberIndices = cell2mat(ZoneColFiberIndices([CellZone(i) AdjacentZones{CellZone(i)}]));

                                ColFiberDistances(AdjacentColFiberIndices) = sqrt((ColFiberPosX(AdjacentColFiberIndices) - CellPosX(i)).^2 + (ColFiberPosY(AdjacentColFiberIndices) - CellPosY(i)).^2);

                                NearbyColFibers = find(ColFiberDistances <= CellRadius);

                                if sum(sum(ColFiberDistribution(NearbyColFibers,:)+FibrinFiberDistribution(NearbyColFibers,:))) > 0

                                    if numel(NearbyColFibers) == 1
                                        nFiberDistribution = (ColFiberDistribution(NearbyColFibers,:)+FibrinFiberDistribution(NearbyColFibers,:))';
                                    else
                                        nFiberDistribution = (sum(ColFiberDistribution(NearbyColFibers,:)+FibrinFiberDistribution(NearbyColFibers,:), 1))';
                                    end

                                    MeannFiberDirX = sum(nFiberDistribution.*cos(2.*FiberDirBins.*pi./180))./sum(nFiberDistribution);
                                    MeannFiberDirY = sum(nFiberDistribution.*sin(2.*FiberDirBins.*pi./180))./sum(nFiberDistribution);
                                    MeannFiberAngle = 1/2*atan2(MeannFiberDirY, MeannFiberDirX);

                                    CFX = cos(MeannFiberAngle);  % collagen fiber cue unit vector x-component
                                    CFY = sin(MeannFiberAngle);  % collagen fiber cue unit vector y-component
                                    CFM = sqrt(MeannFiberDirX^2 + MeannFiberDirY^2);  % collagen fiber cue vector magnitude

                                else

                                    CFX = 0;
                                    CFY = 0;
                                    CFM = 0;

                                end

    %                             clear nFiberDistribution MeannFiberDirX MeannFiberDirY MeannFiberAngle;


                                % compute strength of guidance cues
                                % (normalize and weight the cue
                                % vector magnitudes)
                                Ss = Ws*SM/Ms;
                                Scg = Wcg*CGM/Mcg;
                                Scf = Wcf*CFM/Mcf;
                                Sc = Wc*CM/Mc;


                                % weight the guidance cue unit vectors according to their strength
                                SX = Ss*SX;
                                SY = Ss*SY;
                                CGX = Scg*CGX;
                                CGY = Scg*CGY;
                                CFX = Scf*CFX;
                                CFY = Scf*CFY;
                                CX = Sc*CX;
                                CY = Sc*CY;


                                % compute the resultant guidance cue.  since
                                % stress, fiber, and neighboring cell signals are
                                % bidirectional, need to consider the two
                                % anti-parallel directions.
                                RX = zeros([1 4]);
                                RY = zeros([1 4]);
                                RM = zeros([1 4]);

                                RX(1) = CGX + CX + CFX + SX;
                                RY(1) = CGY + CY + CFY + SY;
                                RM(1) = sqrt(RX(1)^2 + RY(1)^2);
                                RX(2) = CGX + CX - CFX + SX;
                                RY(2) = CGY + CY - CFY + SY;
                                RM(2) = sqrt(RX(2)^2 + RY(2)^2);
                                RX(3) = CGX + CX + CFX - SX;
                                RY(3) = CGY + CY + CFY - SY;
                                RM(3) = sqrt(RX(3)^2 + RY(3)^2);
                                RX(4) = CGX + CX - CFX - SX;
                                RY(4) = CGY + CY - CFY - SY;
                                RM(4) = sqrt(RX(4)^2 + RY(4)^2);


                                % final resultant is whichever of the above cases is the maximum.  if
                                % there is more than one maximum, choose one randomly.
                                k = find(RM == max(RM));
                                k = k(randi([1, numel(k)], 1));


                                % store the resultant angle "THETA" and normalized resultant magnitude "RHO"
                                THETA = atan2(RY(k), RX(k));
                                RHO = RM(k)/(RhoTuningFactor + Scf + Ss + Scg + Sc);
                                SIGMA2 = -2*log(RHO);


                                % compute the probability distribution function for
                                % THETA and RHO and randomly select a direction
                                % from the distribution
                                if RHO >= 0.999

                                    AngleSelection = Directions(Directions == ceil(THETA*180/pi));

                                elseif RHO <= 0.001

                                    AngleSelection = DirBins360(randi([1, numel(DirBins360)], 1));

                                else

                                    WNDistribution = exp(-1.*(DirectionLine*pi/180 - THETA).^2./2./SIGMA2);
                                    for p = 1:1:10
                                        WNDistribution = WNDistribution + exp(-1.*(DirectionLine*pi/180 - THETA + 2*pi*p).^2./2./SIGMA2) + exp(-1.*(DirectionLine*pi/180 - THETA - 2*pi*p).^2./2./SIGMA2);
                                    end
                                    WNDistribution = WNDistribution/sqrt(2*pi*SIGMA2);

                                    WNCumProb = cumtrapz(DirectionLine*pi/180, WNDistribution);
                                    WNBinProb = diff(WNCumProb(1:5:end));
                                    WNBinProb100 = round(WNBinProb/max(WNBinProb)*100);

    %                                         figure(1)
    %                                         clf;
    %                                         plot(DirBins360, WNBinProb100);
    %                                         axis([-180 180 0 100]);
    %                                         pause;

                                    AngleSet = zeros([1 sum(WNBinProb100)]);
                                    counter = 0;
                                    for q = 1:1:numel(DirBins360)
                                        if WNBinProb100(q) == 1
                                            AngleSet(counter + 1) = DirBins360(q);
                                            counter = counter + 1;
                                        elseif WNBinProb100(q) > 1
                                            AngleSet(counter + 1:counter + WNBinProb100(q)) = DirBins360(q);
                                            counter = counter + WNBinProb100(q);
                                        end
                                    end

                                    if counter > 0
                                        AngleSet = AngleSet(randperm(numel(AngleSet)));
                                        AngleSelection = AngleSet(randi([1, numel(AngleSet)], 1));
                                    else
                                        AngleSelection = DirBins360(randi([1, numel(DirBins360)], 1));
                                    end

    %                                 clear WNDistribution WNCumProb WNBinProb WNBinProb100 AngleSet counter p q;

                                end


                                % cell reorientation
                                CellDirX(i) = cos(AngleSelection*pi/180);
                                CellDirY(i) = sin(AngleSelection*pi/180);


                                % cell migration.  cell tries to move a distance=CellSpeed*TimeStep, but
                                % will stop if another cell is in the way
                                CellDistances = sqrt((CellPosX - CellPosX(i)).^2 + (CellPosY - CellPosY(i)).^2);
                                NearbyCells = find(CellDistances > 0 & CellDistances < CellSpeed(i)*TimeStep+2*CellRadius);
                                NearbyCellDir = round(180/pi*atan2(CellPosY(NearbyCells) - CellPosY(i), CellPosX(NearbyCells) - CellPosX(i)));
                                ForwardCells = find(cos((NearbyCellDir - AngleSelection)*pi/180) > 0);

                                if ~isempty(NearbyCells)

                                    J = ceil(5/4*CellSpeed(i)/CellRadius*TimeStep);
                                    stop = 0;
                                    j = 0;

                                    while stop == 0 && j < J

                                        j = j + 1;

                                        Xnew = CellPosX(i) + j/J*CellSpeed(i)*TimeStep*CellDirX(i);
                                        Ynew = CellPosY(i) + j/J*CellSpeed(i)*TimeStep*CellDirY(i);

                                        DistanceNewToNearby = sqrt((CellPosX(NearbyCells(ForwardCells)) - Xnew).^2 + (CellPosY(NearbyCells(ForwardCells)) - Ynew).^2);

                                        if min(DistanceNewToNearby) < 2*CellRadius

                                            j = j - 1;
                                            stop = 1;

                                        end

                                    end
                                    
                                    if j < J
                                    
                                        CellPosX(i) = CellPosX(i) + j/J*CellSpeed(i)*TimeStep*CellDirX(i);
                                        CellPosY(i) = CellPosY(i) + j/J*CellSpeed(i)*TimeStep*CellDirY(i);
                                        CellSpeed(i) = 0;
                                    
                                    else
                                        
                                        CellPosX(i) = CellPosX(i) + CellSpeed(i)*TimeStep*CellDirX(i);
                                        CellPosY(i) = CellPosY(i) + CellSpeed(i)*TimeStep*CellDirY(i);
                                        
                                    end

    %                                 clear J stop j Xnew Ynew DistanceNewToNearby;

                                else

                                    CellPosX(i) = CellPosX(i) + CellSpeed(i)*TimeStep*CellDirX(i);
                                    CellPosY(i) = CellPosY(i) + CellSpeed(i)*TimeStep*CellDirY(i);

                                end

    %                             clear NearbyCells NearbyCellDir ForwardCells;


                                % collagen remodeling
                                AngleSelection = 180/pi*atan(sin(AngleSelection*pi/180)/cos(AngleSelection*pi/180));

                                % collagen fiber rotation.  because the rotation rate depends on the
                                % current fiber orientation, for large rotation rate coefficients,
                                % we break this process into smaller time steps as determined by Q
                                if sum(sum(ColFiberDistribution(NearbyColFibers,:))) > 0

                                    Q = ceil(kColFiberRot(i)*TimeStep/2);

                                    S1indices = find(tan(AngleSelection*pi/180-FiberDirBins*pi/180) > 0 & abs(sin(AngleSelection*pi/180-FiberDirBins*pi/180)) ~= 1);
                                    S2indices = find(tan(AngleSelection*pi/180-FiberDirBins*pi/180) <= 0 & abs(sin(AngleSelection*pi/180-FiberDirBins*pi/180)) ~= 1);
                                    S3indices = find(abs(sin(AngleSelection*pi/180-FiberDirBins*pi/180)) == 1);

                                    NewAngles = zeros(size(FiberDirBins));
                                    NewAngles(S1indices) = FiberDirBins(S1indices) + kColFiberRot(i).*TimeStep./Q.*abs(sin(AngleSelection*pi/180-FiberDirBins(S1indices)*pi/180));
                                    NewAngles(S2indices) = FiberDirBins(S2indices) - kColFiberRot(i).*TimeStep./Q.*abs(sin(AngleSelection*pi/180-FiberDirBins(S2indices)*pi/180));
                                    NewAngles(S3indices) = FiberDirBins(S3indices) + kColFiberRot(i).*TimeStep./Q.*abs(sin(AngleSelection*pi/180-FiberDirBins(S3indices)*pi/180));
                                    NewAngles(NumFiberDirBins+1) = FiberDirBins(S3indices) - kColFiberRot(i).*TimeStep./Q.*abs(sin(AngleSelection*pi/180-FiberDirBins(S3indices)*pi/180));
                                    NewAngles(NewAngles <= -90) = NewAngles(NewAngles <= -90) + 180;
                                    NewAngles(NewAngles > 90) = NewAngles(NewAngles > 90) - 180;

                                    for q = 1:1:Q

                                        NewColFiberDistribution = zeros(size(ColFiberDistribution(NearbyColFibers,:)));

                                        for j = 1:1:NumFiberDirBins+1

                                            if j == S3indices
                                                f = 1/2;
                                            elseif j == NumFiberDirBins+1
                                                f = 1/2;
                                                j = S3indices;
                                            else
                                                f = 1;
                                            end

                                            NewIndex = find(abs(FiberDirBins - NewAngles(j)) == min(abs(FiberDirBins - NewAngles(j))), 1);
                                            OffsetFrac = (NewAngles(j) - FiberDirBins(NewIndex))/5;

                                            if OffsetFrac >= 0

                                                NewColFiberDistribution(:,NewIndex) = NewColFiberDistribution(:,NewIndex) + f*(1-OffsetFrac)*ColFiberDistribution(NearbyColFibers, j);
                                                if NewIndex == NumFiberDirBins
                                                    NewColFiberDistribution(:,1) = NewColFiberDistribution(:,1) + f*OffsetFrac*ColFiberDistribution(NearbyColFibers, j);
                                                else
                                                    NewColFiberDistribution(:,NewIndex+1) = NewColFiberDistribution(:,NewIndex+1) + f*OffsetFrac*ColFiberDistribution(NearbyColFibers, j);
                                                end

                                            else

                                                NewColFiberDistribution(:,NewIndex) = NewColFiberDistribution(:,NewIndex) + f*(1+OffsetFrac)*ColFiberDistribution(NearbyColFibers, j);
                                                if NewIndex == 1
                                                    NewColFiberDistribution(:,NumFiberDirBins) = NewColFiberDistribution(:,NumFiberDirBins) - f*OffsetFrac*ColFiberDistribution(NearbyColFibers, j);
                                                else
                                                    NewColFiberDistribution(:,NewIndex-1) = NewColFiberDistribution(:,NewIndex-1) - f*OffsetFrac*ColFiberDistribution(NearbyColFibers, j);
                                                end

                                            end

                                        end

                                        ColFiberDistribution(NearbyColFibers,:) = NewColFiberDistribution;

                                    end

    %                                 clear j f q Q S1indices S2indices NewAngles NewColFiberDistribution NewIndex OffsetFrac;

                                end

                                % fiber degradation (fibrin degradation only within wound)
                                ColFiberDistribution(NearbyColFibers,:) = ColFiberDistribution(NearbyColFibers,:) - kColFiberDeg(i)*TimeStep*ColFiberDistribution(NearbyColFibers,:);

                                CellPosR = sqrt((CellPosX(i))^2 + (CellPosY(i))^2);
                                CellPosT = atan2(CellPosY(i), CellPosX(i));
                                WoundPerimR = WoundRadiusX*WoundRadiusY/sqrt((WoundRadiusX*sin(CellPosT))^2 + (WoundRadiusY*cos(CellPosT))^2);

                                if WoundPerimR >= CellPosR

                                    FibrinFiberDistribution(NearbyColFibers,:) = FibrinFiberDistribution(NearbyColFibers,:) - kFibFiberDeg(i)*TimeStep*FibrinFiberDistribution(NearbyColFibers,:);

                                end

    %                             clear CellPosR CellPosT WoundPerimR;

                                % collagen fiber deposition
                                if strcmp(deposition, 'AD')

                                    % for aligned deposition
                                    ClosestFiberDirBin = find(abs(FiberDirBins - AngleSelection) == min(abs(FiberDirBins - AngleSelection)), 1);
                                    ColFiberDistribution(NearbyColFibers, ClosestFiberDirBin) = ColFiberDistribution(NearbyColFibers, ClosestFiberDirBin) + kColFiberGen(i)*TimeStep/numel(NearbyColFibers);

%                                             clear ClosestFiberDirBin;

                                else

                                    % for random deposition
                                    RandomFiberDirBin = randi([1 numel(FiberDirBins)], 1);
                                    ColFiberDistribution(NearbyColFibers, RandomFiberDirBin) = ColFiberDistribution(NearbyColFibers, RandomFiberDirBin) + kColFiberGen(i)*TimeStep/numel(NearbyColFibers);

%                                             clear RandomFiberDirBin;

                                end


                                % update collagen fiber patch statistics
                                DirectionBinsTiled = repmat(FiberDirBins', numel(NearbyColFibers), 1);
                                Cbar = sum(ColFiberDistribution(NearbyColFibers,:).*cos(2*DirectionBinsTiled*pi/180), 2)./sum(ColFiberDistribution(NearbyColFibers,:), 2);
                                Sbar = sum(ColFiberDistribution(NearbyColFibers,:).*sin(2*DirectionBinsTiled*pi/180), 2)./sum(ColFiberDistribution(NearbyColFibers,:), 2);

                                ColFiberMA(NearbyColFibers) = 180/pi*1/2*atan2(Sbar,Cbar);
                                ColFiberMVL(NearbyColFibers) = sqrt(Cbar.^2 + Sbar.^2);
                                ColFiberMVX(NearbyColFibers) = ColFiberMVL(NearbyColFibers).*cos(ColFiberMA(NearbyColFibers)*pi/180);
                                ColFiberMVY(NearbyColFibers) = ColFiberMVL(NearbyColFibers).*sin(ColFiberMA(NearbyColFibers)*pi/180);

    %                             clear DirectionBinsTiled Cbar Sbar;

                            end

                        end








                        % eliminate apoptotic cells
                        CellPosX(ApoptosisIndex) = [];
                        CellPosY(ApoptosisIndex) = [];
                        CellDirX(ApoptosisIndex) = [];
                        CellDirY(ApoptosisIndex) = [];
                        CellSpeed(ApoptosisIndex) = [];
                        CellCycleClock(ApoptosisIndex) = [];
                        CellCycleDuration(ApoptosisIndex) = [];
                        ApoptosisClock(ApoptosisIndex) = [];
                        CellZone(ApoptosisIndex) = [];

                        NumCells = NumCells - numel(ApoptosisIndex);

    %                     clear ApoptosisIndex;


                        % update cell zone
                        for i = 1:1:numel(MinX)

                            CellZone(CellPosX >= MinX(i) & CellPosX < MaxX(i) & CellPosY >= MinY(i) & CellPosY < MaxY(i)) = i;

                        end


                        % wrap cells outside tissue boundary
                        CellPosX(CellPosX > xmax) = CellPosX(CellPosX > xmax) - (xmax - xmin);
                        CellPosX(CellPosX < xmin) = CellPosX(CellPosX < xmin) + (xmax - xmin);
                        CellPosY(CellPosY > ymax) = CellPosY(CellPosY > ymax) - (ymax - ymin);
                        CellPosY(CellPosY < ymin) = CellPosY(CellPosY < ymin) + (ymax - ymin);


                        % store data just as the initial state was stored.
                        if mod(t-1,SamplingInterval/TimeStep) == 0

                            ntp = ntp + 1;

                            T(ntp) = (t-1)*TimeStep;

                            NC(ntp) = NumCells;

                            NCF(ntp) = sum(sum(ColFiberDistribution));

                            CPX{ntp} = CellPosX;
                            CPY{ntp} = CellPosY;
                            CDX{ntp} = CellDirX;
                            CDY{ntp} = CellDirY;

                            CFMAtemp = reshape(ColFiberMA, R, C);
                            CFMVLtemp = reshape(ColFiberMVL, R, C);
                            CFMVXtemp = reshape(ColFiberMVX, R, C);
                            CFMVYtemp = reshape(ColFiberMVY, R, C);
                            CFNDtemp = reshape(sum(ColFiberDistribution,2), R, C);

                            FFNDtemp = reshape(sum(FibrinFiberDistribution,2), R, C);

                            for i = 1:AveSpan:R-rem(R,AveSpan)
                                for j = 1:AveSpan:C-rem(C,AveSpan)

                                    Ctemp = sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1).*CFMVLtemp(i:i+AveSpan-1,j:j+AveSpan-1).*cos(2*CFMAtemp(i:i+AveSpan-1,j:j+AveSpan-1)*pi/180)))/sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1)));
                                    Stemp = sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1).*CFMVLtemp(i:i+AveSpan-1,j:j+AveSpan-1).*sin(2*CFMAtemp(i:i+AveSpan-1,j:j+AveSpan-1)*pi/180)))/sum(sum(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1)));
                                    Atemp = 1/2*atan2(Stemp, Ctemp);
                                    Mtemp = sqrt(Stemp^2 + Ctemp^2);

                                    CFMVX((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = Mtemp*cos(Atemp);
                                    CFMVY((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = Mtemp*sin(Atemp);
                                    CFND((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = mean2(CFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1));

                                    FFND((i+AveSpan-1)/AveSpan,(j+AveSpan-1)/AveSpan,ntp) = mean2(FFNDtemp(i:i+AveSpan-1,j:j+AveSpan-1));

    %                                 clear Ctemp Stemp Atemp Mtemp;

                                end
                            end

    %                         clear CFMVXtemp CFMVYtemp CFNDtemp FFNDtemp CFMAtemp CFMVLtemp;


                            CellPosR = sqrt(CellPosX.^2 + CellPosY.^2);
                            CellPosT = atan2(CellPosY, CellPosX);
                            WoundPerimR = WoundRadiusX*WoundRadiusY./sqrt((WoundRadiusX*sin(CellPosT)).^2 + (WoundRadiusY*cos(CellPosT)).^2);

                            WoundCellIndices = find(CellPosR <= WoundPerimR);

                            CWI{ntp,1} = WoundCellIndices;

                            if ~isempty(WoundCellIndices)

                                NCWound(ntp,1) = numel(WoundCellIndices);

                                WoundCellDirX = CellDirX(WoundCellIndices);
                                WoundCellDirY = CellDirY(WoundCellIndices);
                                WoundCellAngle = atan2(WoundCellDirY, WoundCellDirX);

                                MeanWoundCellDirX = mean(cos(2*WoundCellAngle));
                                MeanWoundCellDirY = mean(sin(2*WoundCellAngle));

                                MeanWoundCellAngle(ntp,1) = 180/pi*1/2*atan2(MeanWoundCellDirY, MeanWoundCellDirX);
                                MeanWoundCellLength(ntp,1) = sqrt(MeanWoundCellDirX^2 + MeanWoundCellDirY^2);

                                WoundCellAngle(WoundCellAngle > pi/2 & WoundCellAngle <= pi) = WoundCellAngle(WoundCellAngle > pi/2 & WoundCellAngle <= pi) - pi;
                                WoundCellAngle(WoundCellAngle > -pi & WoundCellAngle <= -pi/2) = WoundCellAngle(WoundCellAngle > -pi & WoundCellAngle <= -pi/2) + pi;
                                WoundCellAngleDistribution(:,ntp,1) = hist(WoundCellAngle*180/pi, FiberDirBins);

                            else

                                NCWound(ntp,1) = 0;
                                MeanWoundCellAngle(ntp,1) = 0;
                                MeanWoundCellLength(ntp,1) = 0;
                                WoundCellAngleDistribution(:,ntp,1) = 0;

                            end

    %                         clear CellPosR CellPosT WoundPerimR WoundCellDirX WoundCellDirY WoundCellAngle MeanWoundCellDirX MeanWoundCellDirY;


                            NCFWound(ntp,1) = sum(sum(ColFiberDistribution(WoundColFiberIndices,:)));
                            NFFWound(ntp,1) = sum(sum(FibrinFiberDistribution(WoundColFiberIndices,:)));

                            if NCFWound(ntp,1) > 0

                                WoundColFiberAngleDistribution(:,ntp,1) = (sum(ColFiberDistribution(WoundColFiberIndices,:), 1))';

                                MeanWoundColFiberDirX = sum(WoundColFiberAngleDistribution(:,ntp,1).*cos(2.*FiberDirBins.*pi./180))./sum(WoundColFiberAngleDistribution(:,ntp,1));
                                MeanWoundColFiberDirY = sum(WoundColFiberAngleDistribution(:,ntp,1).*sin(2.*FiberDirBins.*pi./180))./sum(WoundColFiberAngleDistribution(:,ntp,1));

                                MeanWoundColFiberAngle(ntp,1) = 180/pi*1/2*atan2(MeanWoundColFiberDirY, MeanWoundColFiberDirX);
                                MeanWoundColFiberLength(ntp,1) = sqrt(MeanWoundColFiberDirX^2 + MeanWoundColFiberDirY^2);

                            else

                                MeanWoundColFiberAngle(ntp,1) = 0;
                                MeanWoundColFiberLength(ntp,1) = 0;
                                WoundColFiberAngleDistribution(:,ntp,1) = 0;

                            end

    %                         clear MeanWoundColFiberDirX MeanWoundColFiberDirY;

                        end

                        if mod(t-1,SamplingInterval2/TimeStep) == 0

                            ntp2 = ntp2 + 1;

                            T2(ntp2) = (t-1)*TimeStep;

                            ColPX(:,ntp2) = ColFiberPosX;
                            ColPY(:,ntp2) = ColFiberPosY;
                            ColDX(:,ntp2) = ColFiberMVX;
                            ColDY(:,ntp2) = ColFiberMVY;

                        end

                        % plot
                        %{
                        figure(1)
                        clf;
                        hold on;
                        quiver(ColFiberPosX, ColFiberPosY, ColFiberMVX, ColFiberMVY, 1, 'k', 'LineWidth', 1);
                        quiver(CellPosX, CellPosY, CellDirX, CellDirY, 0.3, 'r', 'LineWidth', 3);
                        plot(CellPosX, CellPosY, 'or', 'MarkerSize', 10, 'LineWidth', 2);
                        plot(WoundPerimX, WoundPerimY, '-b', 'LineWidth', 4);
                        axis(axislimits);
                        axis square;
                        grid on;
                        title('Fibers and Cells', 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                        ylabel('Y Position (um)', 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                        xlabel('X Position (um)', 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                        set(gca, 'LineWidth', 2, 'FontSize', 14, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
                        hold off;
                        %}

    %                     pause;

    %                     display(t);

    %                     display(NumCells);

    %                     toc;

                    end
                    % end of simulation







                    % store results in z variables so they are
                    % saved outside of parfor domain
                    zT{repeat} = T;
                    zNC{repeat} = NC;
                    zNCF{repeat} = NCF;
                    zCWI{repeat} = CWI;
                    zCPX{repeat} = CPX;
                    zCPY{repeat} = CPY;
                    zCDX{repeat} = CDX;
                    zCDY{repeat} = CDY;
                    zCFWI{repeat} = CFWI;
                    zCFMVX{repeat} = CFMVX;
                    zCFMVY{repeat} = CFMVY;
                    zCFND{repeat} = CFND;
                    zFFND{repeat} = FFND;
                    zNCWound{repeat} = NCWound;
                    zMeanWoundCellAngle{repeat} = MeanWoundCellAngle;
                    zMeanWoundCellLength{repeat} = MeanWoundCellLength;
                    zWoundCellAngleDistribution{repeat} = WoundCellAngleDistribution;
                    zNCFWound{repeat} = NCFWound;
                    zMeanWoundColFiberAngle{repeat} = MeanWoundColFiberAngle;
                    zMeanWoundColFiberLength{repeat} = MeanWoundColFiberLength;
                    zWoundColFiberAngleDistribution{repeat} = WoundColFiberAngleDistribution;
                    zNFFWound{repeat} = NFFWound;
                    zT2{repeat} = T2;
                    zColPX{repeat} = ColPX;
                    zColPY{repeat} = ColPY;
                    zColDX{repeat} = ColDX;
                    zColDY{repeat} = ColDY;

                end
%                 matlabpool close;

                save([savename '.mat'], '-mat');

            end

        end

    end

end


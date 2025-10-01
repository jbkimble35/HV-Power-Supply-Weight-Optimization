% This code is to design an inductor based on off-the-shelf cores
% for certain sets of electrical requirements.
%
% Assumptions:
%   - Sine wave on the input and output
%   - Model core loss using standard Steinmetz equation
%   - Model copper loss with full Dowell equations

function y = Ecore_actual_EEER_inductor_LCC(raw,raw1,raw2,raw3,raw4,raw5, ...
    Vin_range, G_range, Po_range, fs_range, Ls_range, Imax_range, Winding_Pattern, ...
    LCC_Q, LCC_f0, LCC_A, LCC_K, LCC_RT, LCC_Ls, LCC_Cs, LCC_Cp, LCC_GT)

% Lowest allowed transformer efficiency
etaInductor = 0.98;

% Max allowable temperature (C)
Tmax = 90;

% Min allowable temperature (C)
Tmin = 25;

% Max allowable current density in the wire (A/m^2)
Jwmax = 500 * 100 * 100;

% Minimal litz diameter one can get (m)
MinLitzDia = 0.05024 / 1000; % AWG44, 0.0316 is AWG48, %0.03983 is AWG46

% Dielectric strength of the insulation material (V/m), discount 50%
dielectricstrength_insulation = 0.5 * 200 * 1000 * 100; % TEFLON

% Minimal air gap (m)
mingap = 10e-6;

% Minimum allowable core cross section radius (m)
MinWinding = 1;

% Maximum turns
MaxWinding = 100;

% Incremental winding
IncreN = 1;

% Maximum layer of winding
MaxMl = 5;

% Incremental layers
IncreMl = 1;

% Minimal wire diameter (m)
MinWireSize = 0.079/1000; % AWG28, 0.35 mm is AWG29, 0.079 is AWG40

% Maximum allowable weight (g)
MaxWeight = 1000;

% g/m^3, density of the core
CoreDensity = 4.8*1000*1000;

% g/m^3, density of copper
CopperDensity = 8.96*1000*1000;

% g/m^3, density of core insulation materials
CoreInsulationDensity = 2.2*1000*1000; % TEFLON

% g/m^3, density of wire insulation materials
WireInsulationDensity = 2.2*1000*1000; % TEFLON

% All discount factors
% Bmax discount factor
BSAT_discount = 0.75;

% Actual core loss is always higher than the calculated
CoreLossMultiple = 1.0;

% Maximum packing factor (copper area compared with total window area)
maxpackingfactor = 0.7;

% Minimum packing factor
minpackingfactor = 0.1;

% Winding factor of litz wire, assuming only 80% of wire size is copper
LitzFactor = 0.8;

% Weight of bobbin as a fraction of the core insulation
BobbinWeightFactor = 0.5;

% Save design results
Design_inductor = zeros(1,42);

% Electrical constants. Normally there is no need to change
% ohmm, resistivity of copper at 100C
rou = 2.34e-8;

% (ohmm), conductivity of copper
sigma = 1/rou;

% H/A·m^2, permeability of free space
u0 = 4 * pi * 10^(-7);

% F/m, permittivity of free space
ebs10 = 8.854e-12;

% MAIN BODY OF THE CODE STARTS FROM HERE

% Read the material properties to get the P at different B and F map from
% the CoreLossData.xlsx file
% This map will be used to calculate core loss later on
[m1, n1] = size(raw1);
XCoreMAT = raw2(2:m1,2);
XCoreFreq = cell2mat(raw1(2:m1,3:n1));
[m1, n1] = size(raw2);
XCoreBfield = cell2mat(raw2(2:m1,3:n1));
[m1, n1] = size(raw3);
XCorePloss = cell2mat(raw3(2:m1,3:n1));
[m1, n1] = size(raw4);
XCoreBSAT = cell2mat(raw4(2:m1,3:n1));
[m1, n1] = size(raw5);
XCoreMU = cell2mat(raw5(2:m1,3:n1));

% Constant
Pbar = 500; % 500 mW/cm3
PFfactor = 1;

% Draw out Pv plot vs B then interpolate
NoMat = m1-1;
Ball = 0.001:0.001:1;
MATcolorvector = rand(NoMat,3);
FreqFlag = zeros(size(1:1:NoMat));

for i = 1:NoMat
    DataSheetFreq = XCoreFreq(i,~isnan(XCoreFreq(i,:)));
    NoFreq = length(DataSheetFreq)/2;
    colorvector = rand(NoFreq,3);
    for j = 1:NoFreq
        % Equation: ConstantA * Bfield + ConstantB
        ConstantA(i,j) = (log10(XCorePloss(i,2*j)) - log10(XCorePloss(i,2*j-1))) / ...
                         (log10(XCoreBfield(i,2*j)) - log10(XCoreBfield(i,2*j-1)));
        ConstantB(i,j) = log10(XCorePloss(i,2*j)) - ConstantA(i,j) * log10(XCoreBfield(i,2*j));
        
        f_at_Pv_500(i,j) = 10^((log10(Pbar) - ConstantB(i,j)) / ConstantA(i,j)); % in T
        f_at_Pv_500(i,j) = DataSheetFreq(2*j-1); % in Hz
        f_at_Pv_500(i,j) = f_at_Pv_500(i,j) * PFfactor;
        
        if (abs(fs_range - f_at_Pv_500(i,j)) / fs_range <= 0.4)
            FreqFlag(i) = 1;
        end
        % Steinmetz
        if (j > 1)
            beta_range(i,j) = log10(XCorePloss(i,2*j) / XCorePloss(i,2*j-1)) / ...
                              log10(XCoreBfield(i,2*j) / XCoreBfield(i,2*j-1));
            % Third point
            XCorePloss_3rd(i,j) = 10^(ConstantA(i,j-1) * log10(XCoreBfield(i,2*j)) + ConstantB(i,j-1));
            alpha_range(i,j) = log10(XCorePloss_3rd(i,j) / XCorePloss(i,2*j)) / ...
                               log10(DataSheetFreq(2*j-3) / DataSheetFreq(2*j-1)); % (f2/f1)^alpha = P2/P1
            K1_range(i,j) = XCorePloss(i,2*j) / (XCoreBfield(i,2*j)^beta_range(i,j) / ...
                           DataSheetFreq(2*j-1)^alpha_range(i,j)); % mW/cm3
            % Repeat frequency 2's steinmetz parameter for frequency 1
            if (j == 2)
                beta_range(i,j-1) = beta_range(i,j);
                alpha_range(i,j-1) = alpha_range(i,j);
                K1_range(i,j-1) = XCorePloss(i,2*j-2) / (XCoreBfield(i,2*j-2)^ ...
                                  beta_range(i,j-1) / DataSheetFreq(2*j-3)^alpha_range(i,j-1));
            end
        end
    end
end

% Core size
[m1, n1] = size(raw);
TransformerCoreIndex = cell2mat(raw(3:m1,1));

XcoreVe = cell2mat(raw(3:m1,3))/(1000^3); % in m
XcoreLe = cell2mat(raw(3:m1,4))/(1000^2);
XcoreAc = cell2mat(raw(3:m1,5))/1000;
XcoreCoreShapeIndex = cell2mat(raw(3:m1,6));
XcorePriW = cell2mat(raw(3:m1,8))/1000;
XcorePriH = cell2mat(raw(3:m1,9))/1000;
XcoreWindowW = cell2mat(raw(3:m1,12))/1000;
XcoreWindowH = 2*cell2mat(raw(3:m1,13))/1000;
ShuffleIndex = 1:1:length(TransformerCoreIndex);

% DESIGN STARTS FROM HERE
CoreMatIndexSweep = find(FreqFlag); % Limit to core materials based on frequency

% Vectorize the design space
[Po, fs, Vin, G, Ls, Imax, matno_record, ShuffleXcoreIndex, Np, Mlp] = ndgrid( ...
    Po_range, fs_range, Vin_range, ...
    G_range, Ls_range, Imax_range, CoreMatIndexSweep, ShuffleIndex, ...
    MinWinding:IncreN:MaxWinding, 1:IncreMl:MaxMl);

Po = reshape(Po,[],1);
fs = reshape(fs,[],1);
Vin = reshape(Vin,[],1);
G = reshape(G,[],1);
Ls = reshape(Ls,[],1);
Imax = reshape(Imax,[],1);
Vinsulation_max = Vin .* G;
matno_record = reshape(matno_record,[],1);
ShuffleXcoreIndex = reshape(ShuffleXcoreIndex,[],1);
Np = reshape(Np,[],1);
Mlp = reshape(Mlp,[],1);
ui = XCoreMU(matno_record);
BSAT = XCoreBSAT(matno_record);

Ve = XcoreVe(ShuffleXcoreIndex); % in cm
Ac = XcoreAe(ShuffleXcoreIndex);
W = XcoreWindowW(ShuffleXcoreIndex);
H = XcoreWindowH(ShuffleXcoreIndex);
Le = XcoreLe(ShuffleXcoreIndex);
PriW = XcorePriW(ShuffleXcoreIndex);
PriH = XcorePriH(ShuffleXcoreIndex);
XcoreIndex = TransformerCoreIndex(ShuffleXcoreIndex);
XcoreCoreShapeIndex = XcoreCoreShapeIndex(ShuffleXcoreIndex);

% Inductance and air gap
Vpri = Vin .* G;
L = Ls;
airgap = u0.*Ac.*Np.^2./L-Le./ui;

% Eliminate some elements based on dimension rule and BSAT rule
KeepAirGap = intersect(find(airgap >= mingap), find(airgap <= 0.2 * Le));
ui = ui ./ (1 + ui .* airgap ./ Le);
Bm_dummy = u0 .* Np .* Imax ./ Le; % T
Keep_Bmindex = find(Bm_dummy < BSAT * BSAT_discount);
KeepIndex = intersect(KeepAirGap, Keep_Bmindex);
Po = Po(KeepIndex);
fs = fs(KeepIndex);
Vin = Vin(KeepIndex);
Vpri = Vpri(KeepIndex);
Imax = Imax(KeepIndex);
Vinsulation_max = Vinsulation_max(KeepIndex);
matno_record = matno_record(KeepIndex);
ui = ui(KeepIndex);
BSAT = BSAT(KeepIndex);

PriW = PriW(KeepIndex);
PriH = PriH(KeepIndex);
H = H(KeepIndex);
W = W(KeepIndex);
Ac = Ac(KeepIndex);
Le = Le(KeepIndex);
Ve = Ve(KeepIndex);
Np = Np(KeepIndex);
Mlp = Mlp(KeepIndex);
L = L(KeepIndex);
airgap = airgap(KeepIndex);
XcoreIndex = XcoreIndex(KeepIndex);
XcoreCoreShapeIndex = XcoreCoreShapeIndex(KeepIndex);

% Find core loss property that's non-zero around the required frequency for each design group
FsnoNonzero = f_at_Pv_500(matno_record,:) > 0;
FsnoIndex = abs(fs - f_at_Pv_500(matno_record,:)) ./ fs <= 0.4;
matfsIndex = FsnoNonzero .* FsnoIndex;
matfs = f_at_Pv_500(matno_record,:) .* matfsIndex;
K1 = K1_range(matno_record,:) .* matfsIndex * 1000; % convert from mW/cm3 to W/m3
alpha = alpha_range(matno_record,:) .* matfsIndex;
beta = beta_range(matno_record,:) .* matfsIndex;
[rowIdcs, colIdcs] = find(matfs > 0);

% Find the indices of unique values in rowIdcs
[UniqueRowIdcs, ind] = unique(rowIdcs, 'rows');
ColDuplicate = sum(matfs(UniqueRowIdcs,:) ~= 0, 2);

% Repeat by the number of loss data of each design point
Po = repelem(Po(UniqueRowIdcs), ColDuplicate);
fs = repelem(fs(UniqueRowIdcs), ColDuplicate);
Vin = repelem(Vin(UniqueRowIdcs), ColDuplicate);
Vpri = repelem(Vpri(UniqueRowIdcs), ColDuplicate);
Imax = repelem(Imax(UniqueRowIdcs), ColDuplicate);
Vinsulation_max = repelem(Vinsulation_max(UniqueRowIdcs), ColDuplicate);
matno_record = repelem(matno_record(UniqueRowIdcs), ColDuplicate);
ui = repelem(ui(UniqueRowIdcs), ColDuplicate);
BSAT = repelem(BSAT(UniqueRowIdcs), ColDuplicate);
PriW = repelem(PriW(UniqueRowIdcs), ColDuplicate);
PriH = repelem(PriH(UniqueRowIdcs), ColDuplicate);
H = repelem(H(UniqueRowIdcs), ColDuplicate);
W = repelem(W(UniqueRowIdcs), ColDuplicate);
Ac = repelem(Ac(UniqueRowIdcs), ColDuplicate);
Le = repelem(Le(UniqueRowIdcs), ColDuplicate);
Ve = repelem(Ve(UniqueRowIdcs), ColDuplicate);
Np = repelem(Np(UniqueRowIdcs), ColDuplicate);
Mlp = repelem(Mlp(UniqueRowIdcs), ColDuplicate);
L = repelem(L(UniqueRowIdcs), ColDuplicate);
airgap = repelem(airgap(UniqueRowIdcs), ColDuplicate);
XcoreIndex = repelem(XcoreIndex(UniqueRowIdcs), ColDuplicate);
XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex(UniqueRowIdcs), ColDuplicate);

% Reformat loss data into one non-zero vector
matfs = nonzeros(reshape(matfs(UniqueRowIdcs,:),[],1));
K1 = nonzeros(reshape(K1(UniqueRowIdcs,:),[],1));
beta = nonzeros(reshape(beta(UniqueRowIdcs,:),[],1));
alpha = nonzeros(reshape(alpha(UniqueRowIdcs,:),[],1));

size(Po)
if (isempty(Po))
    y = 0;
else
    % Repeat elements by Primary Wire Number of Strands
    skindepth = 1 ./ sqrt(pi * fs .* u0 / rou);
    ds = max(skindepth, MinLitzDia * ones(size(skindepth))); % take the skin depth litz
    ds = MinLitzDia * ones(size(skindepth));
    MinPriNstrands = floor((Imax ./ Jwmax) ./ (pi * ds.^2 / 4)) + 1;
    MaxPriNstrands = floor((Imax ./ Jwmax*1.0) ./ (pi * ds.^2 / 4)) + 1;

    Po = repelem(Po, [MaxPriNstrands - MinPriNstrands + 1]);
    fs = repelem(fs, [MaxPriNstrands - MinPriNstrands + 1]);
    Vin = repelem(Vin, [MaxPriNstrands - MinPriNstrands + 1]);
    Vpri = repelem(Vpri, [MaxPriNstrands - MinPriNstrands + 1]);
    Imax = repelem(Imax, [MaxPriNstrands - MinPriNstrands + 1]);
    Vinsulation_max = repelem(Vinsulation_max, [MaxPriNstrands - MinPriNstrands + 1]);
    matno_record = repelem(matno_record, [MaxPriNstrands - MinPriNstrands + 1]);
    ui = repelem(ui, [MaxPriNstrands - MinPriNstrands + 1]);
    BSAT = repelem(BSAT, [MaxPriNstrands - MinPriNstrands + 1]);
    PriW = repelem(PriW, [MaxPriNstrands - MinPriNstrands + 1]);
    PriH = repelem(PriH, [MaxPriNstrands - MinPriNstrands + 1]);
    H = repelem(H, [MaxPriNstrands - MinPriNstrands + 1]);
    W = repelem(W, [MaxPriNstrands - MinPriNstrands + 1]);
    Ac = repelem(Ac, [MaxPriNstrands - MinPriNstrands + 1]);
    Le = repelem(Le, [MaxPriNstrands - MinPriNstrands + 1]);
    Ve = repelem(Ve, [MaxPriNstrands - MinPriNstrands + 1]);
    Np = repelem(Np, [MaxPriNstrands - MinPriNstrands + 1]);
    Mlp = repelem(Mlp, [MaxPriNstrands - MinPriNstrands + 1]);
    matfs = repelem(matfs, [MaxPriNstrands - MinPriNstrands + 1]);
    K1 = repelem(K1, [MaxPriNstrands - MinPriNstrands + 1]);
    beta = repelem(beta, [MaxPriNstrands - MinPriNstrands + 1]);
    alpha = repelem(alpha, [MaxPriNstrands - MinPriNstrands + 1]);
    L = repelem(L, [MaxPriNstrands - MinPriNstrands + 1]);
    airgap = repelem(airgap, [MaxPriNstrands - MinPriNstrands + 1));
    XcoreIndex = repelem(XcoreIndex, [MaxPriNstrands - MinPriNstrands + 1]);
    XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex, [MaxPriNstrands - MinPriNstrands + 1]);
    Pri_Nstrands = repmat((MinPriNstrands(1):1:MaxPriNstrands(1)), length(MaxPriNstrands), 1);
    
    % Calculate core loss (W)
    Iprms = Imax ./ sqrt(2);
    ue = ui ./ (1 + ui .* airgap ./ Le);
    Bm = u0 .* Np .* Imax ./ Le ./ ue; % T
    % Window area
    Wa = H .* W;
    % Core volume (m3)
    Vcore = Ve;
    % Core weight (g)
    Wcore = Vcore .* CoreDensity;
    % Check core loss
    Pcore = CoreLossMultiple .* Vcore .* K1 .* fs.^alpha .* Bm.^beta;
    
    % Wire
    % Recalculate ds
    skindepth = 1 ./ sqrt(pi * fs * u0 / rou);
    ds = max(skindepth, MinLitzDia * ones(size(skindepth))); % take the skin depth litz
    ds = MinLitzDia * ones(size(skindepth)); % take the smallest litz
    % Primary wire diameter (m)
    Pri_WireSize = sqrt(Pri_Nstrands * pi * ds.^2 / 4 / LitzFactor / pi) * 2;
    % Primary wire diameter (m) including the insulation layer
    Pri_FullWireSize = Pri_WireSize + (Vin ./ dielectricstrength_insulation) * 2;
    % Packing factors
    CopperPacking = (pi * Pri_WireSize.^2 * Np / 4) ./ (H .* W);
    OverallPacking = (pi * Pri_FullWireSize.^2 * Np / 4) ./ (H .* W);

    % Winding structures
    % Core insulation thickness needed
    CoreInsulationThickness = Vinsulation_max ./ dielectricstrength_insulation;
    % Turns per layer
    Pri_PerLayer = floor(Np ./ Mlp);
    % Total length of windings
    % For XcoreCoreShapeIndex == 1, EE cores
    TLP = Np * 2 * (PriW + PriH + 4 * CoreInsulationThickness + 2 * Mlp .* Pri_FullWireSize);
    % Recalculate XcoreCoreShapeIndex == 2, ER cores
    SelecIndex = find(XcoreCoreShapeIndex == 2);
    TLP(SelecIndex) = 2 * pi * Np(SelecIndex) .* (PriW(SelecIndex) / 2 + ...
    CoreInsulationThickness(SelecIndex) + 0.5 * Mlp(SelecIndex) .* Pri_FullWireSize(SelecIndex));
    
    % Copper Loss
    PriKlayer = sqrt(pi * Pri_Nstrands) .* ds / 2 ./ (Pri_WireSize);
    Pri_xp = ds ./ 2 ./ skindepth .* sqrt(pi * PriKlayer);
    Pri_Rdc = rou.*TLP./(pi.*Pri_WireSize.^2./4);
    Pri_Fr = Pri_xp .* ((sin(2 .* Pri_xp) + sin(2 .* Pri_xp)) ./ (cosh(2 .* Pri_xp) - cos(2 .* Pri_xp)) + ...
        2 + (Mlp .* 2 .* Pri_Nstrands - 1) ./ 3 .* (sinh(Pri_xp) - sin(Pri_xp)) ./ ...
        (cosh(Pri_xp) - cos(Pri_xp)));
    Pri_Rac = Pri_Rdc .* Pri_Fr;
    Pcopper = (Iprms.^2) .* Pri_Rac;
    
    % Calculate temperature rise
    Rth = 16.31 * Le.^-3 .* (Ac .* Wa).^(-0.405); % what is Wa??
    Tafterloss = Rth .* (Pcopper + Pcore) + 25;
    
    % Calculate the weight
    WeightPri_copper = pi * Pri_WireSize.^2 / 4 .* TLP .* CopperDensity;
    WeightPri_Insu = pi * (Pri_FullWireSize.^2 - Pri_WireSize.^2) / 4 .* TLP .* WireInsulationDensity;
    % For XcoreCoreShapeIndex == 1, EE cores
    WeightCore_Insu = (2 .* H .* (PriW + 2 * PriH) + 4 .* W .* (PriW + 2 * PriH) + ...
        H .* (2 * PriW + 2 * PriH)) .* CoreInsulationThickness .* CoreInsulationDensity;
    % Recalculate XcoreCoreShapeIndex == 2, ER cores
    SelecIndex = find(XcoreCoreShapeIndex == 2);
    WeightCore_Insu(SelecIndex) = (sqrt(2) * pi * W(SelecIndex) .* PriW(SelecIndex) + ...
        sqrt(2) * pi * 2 * W(SelecIndex) .* PriW(SelecIndex) + ...
        H(SelecIndex) .* pi .* PriW(SelecIndex)) .* CoreInsulationThickness(SelecIndex) .* CoreInsulationDensity;

    TotalWeight = Wcore + WeightPri_copper + WeightCore_Insu;
    
    % Filter the designs
    B_index = find(Bm <= BSAT * BSAT_discount);
    P_loss_index = find((Pcopper + Pcore) <= Po * (1 - etaInductor));
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);
    TotalWeight_index = find(TotalWeight <= MaxWeight);
    
    OverallPackingmax_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmin_index = find(OverallPacking <= maxpackingfactor);
    Mlp_index = find(Mlp .* Pri_FullWireSize <= W - 2 * CoreInsulationThickness);
    PerLayer_index = find(Pri_PerLayer .* Pri_FullWireSize <= H - 2 * CoreInsulationThickness);
    
    Index_Meet_All = intersect(B_index, P_loss_index);
    Index_Meet_All = intersect(Index_Meet_All, Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All, Tmin_index);
    Index_Meet_All = intersect(Index_Meet_All, TotalWeight_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmax_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmin_index);
    Index_Meet_All = intersect(Index_Meet_All, Mlp_index);
    Index_Meet_All = intersect(Index_Meet_All, PerLayer_index);
    % Sort by total weight and keep only the lightest one
    [WeightSort, SortIndex] = sort(TotalWeight(Index_Meet_All));
    if (length(SortIndex) >= 1)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1:1));
    
        Design_inductor(:,1) = Po(TotalWeightSortIndex);
        Design_inductor(:,2) = Vin(TotalWeightSortIndex);
        Design_inductor(:,3) = Vpri(TotalWeightSortIndex);
        
        Design_inductor(:,4) = Vinsulation_max(TotalWeightSortIndex);
        Design_inductor(:,5) = fs(TotalWeightSortIndex);
        Design_inductor(:,6) = matno_record(TotalWeightSortIndex);
        Design_inductor(:,7) = matfs(TotalWeightSortIndex);
        
        Design_inductor(:,8) = PriW(TotalWeightSortIndex);
        Design_inductor(:,9) = PriH(TotalWeightSortIndex);
        Design_inductor(:,10) = Ac(TotalWeightSortIndex);
        Design_inductor(:,11) = H(TotalWeightSortIndex);
        Design_inductor(:,12) = W(TotalWeightSortIndex);
        Design_inductor(:,13) = Np(TotalWeightSortIndex);
        Design_inductor(:,14) = Bm(TotalWeightSortIndex);
        
        Design_inductor(:,15) = Pri_WireSize(TotalWeightSortIndex);
        Design_inductor(:,16) = Pri_FullWireSize(TotalWeightSortIndex);
        Design_inductor(:,17) = Imax(TotalWeightSortIndex) ./ ...
            (pi * Pri_Nstrands(TotalWeightSortIndex) .* ds(TotalWeightSortIndex).^2 / 4);
        Design_inductor(:,18) = Pri_Nstrands(TotalWeightSortIndex);
        Design_inductor(:,19) = Pri_PerLayer(TotalWeightSortIndex);
        Design_inductor(:,20) = Mlp(TotalWeightSortIndex);
        
        Design_inductor(:,21) = CopperPacking(TotalWeightSortIndex);
        Design_inductor(:,22) = OverallPacking(TotalWeightSortIndex);
        Design_inductor(:,23) = Pcore(TotalWeightSortIndex);
        Design_inductor(:,24) = Pcopper(TotalWeightSortIndex);
        Design_inductor(:,25) = Wcore(TotalWeightSortIndex);
        Design_inductor(:,26) = WeightPri_copper(TotalWeightSortIndex);
        Design_inductor(:,27) = WeightPri_Insu(TotalWeightSortIndex);
        Design_inductor(:,28) = WeightCore_Insu(TotalWeightSortIndex);
        Design_inductor(:,29) = TotalWeight(TotalWeightSortIndex);
        Design_inductor(:,30) = Tafterloss(TotalWeightSortIndex);
        Design_inductor(:,31) = L(TotalWeightSortIndex);
        Design_inductor(:,32) = airgap(TotalWeightSortIndex);
        Design_inductor(:,33) = XcoreIndex(TotalWeightSortIndex);
        Design_inductor(:,34) = LCC_Q * ones(size(TotalWeightSortIndex));
        Design_inductor(:,35) = LCC_f0 * ones(size(TotalWeightSortIndex));
        Design_inductor(:,36) = LCC_A * ones(size(TotalWeightSortIndex));
        Design_inductor(:,37) = LCC_K * ones(size(TotalWeightSortIndex));
        Design_inductor(:,38) = LCC_RT * ones(size(TotalWeightSortIndex));
        Design_inductor(:,39) = LCC_Ls * ones(size(TotalWeightSortIndex));
        Design_inductor(:,40) = LCC_Cs * ones(size(TotalWeightSortIndex));
        Design_inductor(:,41) = LCC_Cp * ones(size(TotalWeightSortIndex));
        Design_inductor(:,42) = LCC_GT * ones(size(TotalWeightSortIndex));

        y = Design_inductor;
    else
        y = zeros(1,42);
    end
end
% This code is to design an inductor based on off-the-shelf cores
% for certain sets of electrical requirements.
% Assumptions:
%   - Sine wave on the input and output
%   - Model core loss using standard Steinmetz equation
%   - Model copper loss with full Dowell equations

function y = Ecore_actual_EEER_inductor_LCC(raw,raw1,raw2,raw3,raw4,raw5,raw6, ...
    Vin_range, G_range, Po_range, fs_range, Ls_range, Imax_range, Winding_Pattern, ...
    LCC_Q, LCC_f0, LCC_A, LCC_K, LCC_RT, LCC_Ls, LCC_Cs, LCC_Cp, LCC_GT) %#ok<INUSD>

% Winding_Pattern is not used for some reason, but I think it corresponds
% to the XcoreCoreShapeIndex, I might be missing something.

%% Tunable Parameters

% Lowest allowed transformer efficiency
etaInductor = 0.9;
% Max allowable temperature (C)
Tmax = 100;
% Min allowable temperature (C)
Tmin = 25;
% Max allowable current density in the wire (A/m^2)
Jwmax = 500 * 100 * 100;

% Minimal litz diameter one can get (m)
MinLitzDia = 0.05024 / 1000; % AWG44
% 0.0316 is AWG48, 0.03983 is AWG46

% Dielectric strength of the insulation material (V/m), discount 50%
dielectricstrength_insulation = 0.5 * 200 * 1000 * 100; % TEFLON
% Minimal air gap (m)
mingap = 0;
% Minimum allowable core cross section radius (m)
MinWinding = 1;
% Maximum turns
MaxWinding = 2000;
% Incremental winding
IncreN = 1;
% Maximum layer of winding
MaxMl = 20;
% Incremental layers. The layers of a transformer reference each wrap of
% turns that fills the window height before moving on to the next level.
% Once one layer fills, the next layer is wound on top, seperated by an
% insulation layer.
IncreMl = 1;



% Minimal wire diameter (m)
% Not used here for some reason
MinWireSize = 0.079/1000; %#ok<NASGU> % AWG28, 0.35 mm is AWG29, 0.079 is AWG40
% Maximum allowable weight (g)
MaxWeight = 2000;
% g/m^3, density of copper
CopperDensity = 8.96*1000*1000;
% g/m^3, density of core insulation materials
CoreInsulationDensity = 2.2*1000*1000; % TEFLON
% g/m^3, density of wire insulation materials
WireInsulationDensity = 2.2*1000*1000; % TEFLON

% All discount factors
% Bmax discount factor
BSAT_discount = 0.9;
% Actual core loss is always higher than the calculated
CoreLossMultiple = 1.0;
% Maximum packing factor (copper area compared with total window area)
maxpackingfactor = 0.7;
% Minimum packing factor
minpackingfactor = 0.001;
% Winding factor of litz wire, assuming only 80% of wire size is copper
LitzFactor = 0.8;
% Weight of bobbin as a fraction of the core insulation
% Not used here for some reason
BobbinWeightFactor = 0.5; %#ok<NASGU>

%% Electrical Parameters

% Preallocation of design table
Design_inductor = zeros(1,43);
% Electrical constants. Normally there is no need to change
% ohm*m, resistivity of copper at 100C
rou = 2.3*1e-8;
% (ohm*m), conductivity of copper
% Not used here for some reason
sigma = 1/rou; %#ok<NASGU>
% H/A·m^2, permeability of free space
u0 = 4*pi*10^(-7);
% F/m, permittivity of free space
% Not used here for some reason
ebs10 = 8.854*1e-12; %#ok<NASGU>

%%Function Body

% Parse core loss material maps from CoreLossData.xlsx
[m1,n1]   = size(raw1);
XCoreFreq = cell2mat(raw1(2:m1,3:n1));
[m1,n1]   = size(raw2);
XCoreBfield = cell2mat(raw2(2:m1,3:n1));
[m1,n1]   = size(raw3);
XCorePloss   = cell2mat(raw3(2:m1,3:n1));
% Only parses column C
[m1,~]    = size(raw4);
XCoreBSAT = cell2mat(raw4(2:m1,3));
% Only parses column C for now.
[m1,~]    = size(raw5);
XCoreMU   = cell2mat(raw5(2:m1,3));
[m1,~]    = size(raw6);
CoreDensity   = cell2mat(raw6(2:m1,3))*1000000;


% Reference loss level and PF-factor for Steinmetz equations
Pbar = 500;      % mW/cm^3
PFfactor = 1;

% Build Steinmetz parameter sets around the target frequency
% Number of materials excluding the header column
NoMat = m1-1;
% Preallocates flag that will be 0 if the material doesn't have a frequency
% near the design point, and 1 if it does, for each of the material indexes
FreqFlag = zeros(size(1:1:NoMat));
% Sweeps all of the materials. This is the core loss curve fitting loop for
% each material. It processes the datasheet parameters into Steinmetz parameters.
for i = 1:1:NoMat
    % Obtains all frequency values for the material, ignoring NaN values.
    DataSheetFreq = XCoreFreq(i,~isnan(XCoreFreq(i,:)));
    % This defines the number of B-CL flux density vs. loss data sets per
    % frequency, in this case being 2.
    NoFreq = length(DataSheetFreq)/2;

    % Loops over each frequency dataset j for each material i. If desiring
    % the use of more B-CL pairs per frequency, a linear regression fit is
    % needed instead of this 2-point slope fit.

    for j = 1:1:NoFreq
        % log10(P) = A*log10(B) + B0 at this data-sheet freq
        
        % A is the slope, and B is the y-intercept of the line between the
        % 2 B-CL values.
        ConstantA(i,j) = (log10(XCorePloss(i,2*j)) - log10(XCorePloss(i,2*j-1)))/( ...
            log10(XCoreBfield(i,2*j)) - log10(XCoreBfield(i,2*j-1)));
        ConstantB(i,j) = log10(XCorePloss(i,2*j)) - ConstantA(i ,j)*log10(XCoreBfield( ...
            i,2*j));

        % Loss is normalized to the reference power density of 500 mW/cm^3.
        % First it finds flux density where losses equal the reference
        % level, then saves that frequency, then forms stress index
        % combining B and f.
        B_atPv_500(i,j) = 10^((log10(Pbar) - ConstantB(i,j))/ConstantA(i,j)); % in T
        F_atPv_500(i,j) = DataSheetFreq(2*j-1); % in Hz
        PF_atPv_500(i,j) = B_atPv_500(i,j)*F_atPv_500(i,j)^PFfactor;

        % This checks if material data is 0.X away from the desired
        % frequency, and if it is, FreqFlag is raised for that material.
        if (abs(fs_range - F_atPv_500(i,j))/fs_range <= 0.4)
            FreqFlag(i) = 1;
        end

        % Steinmetz
        if (j > 1)

            % Steinmetz exponents across two adjacent frequencies
            beta_range(i,j) = log10(XCorePloss(i,2*j)/XCorePloss(i,2*j-1))/log10( ...
                XCoreBfield(i,2*j)/XCoreBfield(i,2*j-1));
            % Extrapolate third point on previous line (same B2)
            XCorePloss_3rd(i,j) = 10.^(ConstantA(i,j-1)*log10(XCoreBfield(i,2*j))+ ...
                ConstantB(i,j-1));
            % alpha exponent from frequency dependence
            alpha_range(i,j) = log10(XCorePloss_3rd(i,j)/XCorePloss(i,2*j))/log10( ...
                DataSheetFreq(2*j-3)/DataSheetFreq(2*j-1)); %(f2/f1)^alpha = P2/P1;
            % K coefficient in Steinmetz equation
            K1_range(i,j) = XCorePloss(i,2*j)/(XCoreBfield(i,2*j)^beta_range(i,j))/( ...
                DataSheetFreq(2*j-1)^alpha_range(i,j)); %mW/cm3

            % Populate j-1 if missing
            if j == 2
                beta_range(i,j-1)  = beta_range(i,j);
                alpha_range(i,j-1) = alpha_range(i,j);
                K1_range(i,j-1) = XCorePloss(i,2*j-2)/(XCoreBfield(i,2*j-2)^ ...
                    beta_range(i,j-1))/(DataSheetFreq(2*j-3)^alpha_range(i,j-1));            
            end
        end
    end
end


% Core size data parsing
[m1,~] = size(raw);
TransformerCoreIndex = cell2mat(raw(2:m1,1));
% Effective core volume in m^3
XcoreVe = cell2mat(raw(2:m1,3))/(1000^3);
% Cross-sectional area in m^2
XcoreAe = cell2mat(raw(2:m1,4))/(1000^2);
% Magnetic path length in m
XcoreLe = cell2mat(raw(2:m1,5))/1000;

% 1 for EE, 2 for ER, will be 3 for U and 4 for UR
XcoreCoreShapeIndex = cell2mat(raw(2:m1,6));
% Primary winding window dimensions in m
XcorePriW = cell2mat(raw(2:m1,8))/1000;
XcorePriH = cell2mat(raw(2:m1,9))/1000;
% Not used here for some reason
XcoreSecW = cell2mat(raw(2:m1,10))/1000; %#ok<NASGU>
% Not used here for some reason
XcoreSecH = cell2mat(raw(2:m1,11))/1000; %#ok<NASGU>


% PriW means primary winding width
% PriH means primary winding height
% SecW means secondary winding width
% SecH means secondary winding height

% If EE core:
    % If center-leg winding pattern:
        % PriW=Lc and PriH=T
        % SecW=N/A and SecH=N/A
    % If double-leg winding pattern:
        % PriW=Lc and PriH=T
        % SecW=0.5Lc and SecH=T

% If ER core:
    % If center-leg:
        % PriW=rAc and PriH=N/A
        % SecW and SecH N/A
    % If double-leg:
        % I don't think this is in the script

% Lc is center-leg width
% T is core thickness in EE/U cores
% rAc is core radius in ER, UR cores.

XcoreWindowW = cell2mat(raw(2:m1,12))/1000;
XcoreWindowH = cell2mat(raw(2:m1,13))/1000;
% removed 2x for windowH
ShuffleIndex = 1:1:length(TransformerCoreIndex);

%% Design Starts Here

% Material indices with data near fs
CoreMatIndexSweep = find(FreqFlag);  

% Need to make this more efficient somehow. Currently eats up too much RAM
% Creates 10-D array of all ranges of values for sweeping.
[Po, fs, Vin, G, Ls, Imax, matno_record, ShuffleXcoreIndex, Np, Mlp] = ndgrid( ...
    Po_range, fs_range, Vin_range, G_range, Ls_range, Imax_range, ...
    CoreMatIndexSweep, ShuffleIndex, MinWinding:IncreN:MaxWinding, 1:IncreMl:MaxMl);

% Flattens the 10-D array into a very long column vector
Po = reshape(Po,[],1);
fs = reshape(fs,[],1);
Vin = reshape(Vin,[],1);
G = reshape(G,[],1);
Ls = reshape(Ls,[],1);
Imax = reshape(Imax,[],1);
matno_record = reshape(matno_record,[],1);
ShuffleXcoreIndex = reshape(ShuffleXcoreIndex,[],1);
Np = reshape(Np,[],1);
Mlp = reshape(Mlp,[],1);

% Maximum voltage stress for each candidate
Vinsulation_max = Vin.*G;

% Finds BSAT and relative permittivity for each candidate design, using matno_record
% as the index deciding which material is assigned to which index.
ui = XCoreMU(matno_record);
BSAT = XCoreBSAT(matno_record);

% ShuffleXcoreIndex is a very large column vector containing values of indices
% for rows in the CoreSizeData file, and in the following lines, creates very large column vectors
% corresponding to the dimensions for each core geometry chosen. Because
% each core geometry consists of 8 parameters, this is just unfolding those.
W       = XcoreWindowW(ShuffleXcoreIndex);
H       = XcoreWindowH(ShuffleXcoreIndex);
PriW    = XcorePriW(ShuffleXcoreIndex);
PriH    = XcorePriH(ShuffleXcoreIndex);
Ve      = XcoreVe(ShuffleXcoreIndex);
Ac      = XcoreAe(ShuffleXcoreIndex);
Le      = XcoreLe(ShuffleXcoreIndex);
XcoreIndex = TransformerCoreIndex(ShuffleXcoreIndex);
XcoreCoreShapeIndex = XcoreCoreShapeIndex(ShuffleXcoreIndex);


% Primary winding voltage seen by inductor
Vpri    = Vin.*G;
% Shorthand for ease of readability
L       = Ls;
% Effective air gap length based on the magnetic circuit relation between
% inductance, turns, effective length, cross-sectional area, etc.
% SOMETHING IS WRONG WITH THE AIRGAP COMPUTATION, I BELIEVE THE
% CORESIZEDATA FILE IS SOMEHOW BEING READ WRONG OR THE EQUATIONS ARE WRONG.
% MOST LIKELY THE ETD AND ER CORES ARE BEING TREATED LIKE EE CORES AND VICE
% VERSA.

% Safety check for air gaps. If the computed air gap is negative, it means
% the inductance isn't enough.

airgap = u0.*Ac.*Np.^2./L - Le./ui;
% Eliminate unrealistic elements: require a feasible positive gap
KeepAirGap = intersect(find(airgap >= mingap), find(airgap <= 0.2.*Le));


% Optional: check if no candidates survive
if isempty(KeepAirGap)
    warning("No feasible airgap: requested inductance larger than ungapped L0 for all cores/turns.");
end

% Effective permeability of the gapped core
ue = ui./(1+ui.*airgap./Le);
% Estimates peak flux density
Bm_dummy = u0.*Np.*Imax./Le.*ue;
% Keeps designs where saturation > estimate with safety margin
Keep_Bmindex = find(Bm_dummy < BSAT*BSAT_discount);
% Keeps designs that pass both airgap and flux density
KeepIndex = intersect(KeepAirGap, Keep_Bmindex);

% NEW AIRGAP CHECK
L_calc = u0.*Np.^2.*Ac.*ui./(Le + ui.*airgap);
L_index = find(abs(L_calc - L)./L <= 0.05);
KeepIndex = intersect(KeepIndex, L_index);

% Debug
if isempty(KeepIndex)
    warning('No Successful Indexes');
    y = Design_inductor;
    return;
end

% Keeps all resulting indices
Po                 = Po(KeepIndex);
fs                 = fs(KeepIndex);
Vin                = Vin(KeepIndex);
Vpri               = Vpri(KeepIndex);
Imax               = Imax(KeepIndex);
Vinsulation_max    = Vinsulation_max(KeepIndex);
matno_record       = matno_record(KeepIndex);
ui                 = ui(KeepIndex);
BSAT               = BSAT(KeepIndex);
PriW               = PriW(KeepIndex);
PriH               = PriH(KeepIndex);
H                  = H(KeepIndex);
W                  = W(KeepIndex);
Ac                 = Ac(KeepIndex);
Le                 = Le(KeepIndex);
Ve                 = Ve(KeepIndex);
Np                 = Np(KeepIndex);
Mlp                = Mlp(KeepIndex);
L                  = L(KeepIndex);
airgap             = airgap(KeepIndex);
XcoreIndex         = XcoreIndex(KeepIndex);
XcoreCoreShapeIndex= XcoreCoreShapeIndex(KeepIndex);

% Filters out <=0 core loss at all frequencies, material frequencies >20% away from fs,
% and updates candidate frequency points into matfs.
FsnoNonzero = F_atPv_500(matno_record,:) > 0;
FsnoIndex = abs(fs - F_atPv_500(matno_record,:)) ./ fs <= 0.2;
matfsIndex = FsnoNonzero .* FsnoIndex;
matfs = F_atPv_500(matno_record,:) .* matfsIndex;

% Converts K1 from mW/cm3 to W/m3 and builds core-loss steinmetz parameters
% K1, alpha, and beta for the frequencies given by matfs. Then rowIdcs
% stores the indices in the arrays that correspond to valid frequency
% points.
K1 = K1_range(matno_record,:) .* matfsIndex * 1000;
alpha = alpha_range(matno_record,:) .* matfsIndex;
beta = beta_range(matno_record,:) .* matfsIndex;
[rowIdcs, ~] = find(matfs > 0);

% Finds the indices of unique values in rowIdcs, makes col vector where
% each entry is "how many nonzero freq.s exist in a given row of matfs"
% i.e. how many frequencies per each material.
[UniqueRowIdcs, ~] = unique(rowIdcs, 'rows');
ColDuplicate = sum(matfs(UniqueRowIdcs,:) ~= 0, 2);


% Repeat by the number of loss data of each design point
Po                  = repelem(Po(UniqueRowIdcs), ColDuplicate);
fs                  = repelem(fs(UniqueRowIdcs), ColDuplicate);
Vin                 = repelem(Vin(UniqueRowIdcs), ColDuplicate);
Vpri                = repelem(Vpri(UniqueRowIdcs), ColDuplicate);
Imax                = repelem(Imax(UniqueRowIdcs), ColDuplicate);
Vinsulation_max     = repelem(Vinsulation_max(UniqueRowIdcs), ColDuplicate);
matno_record        = repelem(matno_record(UniqueRowIdcs), ColDuplicate);
ui                  = repelem(ui(UniqueRowIdcs), ColDuplicate);
BSAT                = repelem(BSAT(UniqueRowIdcs), ColDuplicate);
PriW                = repelem(PriW(UniqueRowIdcs), ColDuplicate);
PriH                = repelem(PriH(UniqueRowIdcs), ColDuplicate);
H                   = repelem(H(UniqueRowIdcs), ColDuplicate);
W                   = repelem(W(UniqueRowIdcs), ColDuplicate);
Ac                  = repelem(Ac(UniqueRowIdcs), ColDuplicate);
Le                  = repelem(Le(UniqueRowIdcs), ColDuplicate);
Ve                  = repelem(Ve(UniqueRowIdcs), ColDuplicate);
Np                  = repelem(Np(UniqueRowIdcs), ColDuplicate);
Mlp                 = repelem(Mlp(UniqueRowIdcs), ColDuplicate);
L                   = repelem(L(UniqueRowIdcs), ColDuplicate);
airgap              = repelem(airgap(UniqueRowIdcs), ColDuplicate);
XcoreIndex          = repelem(XcoreIndex(UniqueRowIdcs), ColDuplicate);
XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex(UniqueRowIdcs), ColDuplicate);

% Reformat loss data into one non-zero vector
matfs = nonzeros(reshape(matfs(UniqueRowIdcs,:)',[],1));
K1 = nonzeros(reshape(K1(UniqueRowIdcs,:)',[],1));
beta = nonzeros(reshape(beta(UniqueRowIdcs,:)',[],1));
alpha = nonzeros(reshape(alpha(UniqueRowIdcs,:)',[],1));

size(Po)
if (isempty(Po))
    warning('No successful inductor results Test1, Inductor.');
    y = 0;
else
    % Repeat elements by Primary Wire Number of Strands
    skindepth = 1./sqrt(pi*fs*u0/rou);
    % Overwritten variable
    ds = max(skindepth, MinLitzDia*ones(size(skindepth))); % take the skin depth litz

    MinPriNstrands = floor((Imax./Jwmax)./(pi*ds.^2/4))+1;
    MaxPriNstrands = floor((Imax./Jwmax*1.0)./(pi*ds.^2/4))+1;

    % Test
    priCount = MaxPriNstrands - MinPriNstrands + 1;
    if any(priCount < 1)
        bad = find(priCount < 1, 1);
        error('Strands:BadPrimaryCounts', ...
              'Computed primary strand counts <=0 at row %d. Check Po, Vppeak, Jwmax, MinLitzDia.', bad);
    end
    % Testend
        
    Po                  = repelem(Po, priCount);
    fs                  = repelem(fs, priCount);
    Vin                 = repelem(Vin, priCount);
    Vpri                = repelem(Vpri, priCount);
    Imax                = repelem(Imax, priCount);
    Vinsulation_max     = repelem(Vinsulation_max, priCount);
    matno_record        = repelem(matno_record, priCount);
    ui                  = repelem(ui, priCount);
    BSAT                = repelem(BSAT, priCount);
    PriW                = repelem(PriW, priCount);
    PriH                = repelem(PriH, priCount);
    H                   = repelem(H, priCount);
    W                   = repelem(W, priCount);
    Ac                  = repelem(Ac, priCount);
    Le                  = repelem(Le, priCount);
    Ve                  = repelem(Ve, priCount);
    Np                  = repelem(Np, priCount);
    Mlp                 = repelem(Mlp, priCount);
    matfs               = repelem(matfs, priCount);
    K1                  = repelem(K1, priCount);
    beta                = repelem(beta, priCount);
    alpha               = repelem(alpha, priCount);
    L                   = repelem(L, priCount);
    airgap              = repelem(airgap, priCount);
    XcoreIndex          = repelem(XcoreIndex, priCount);
    XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex, priCount);


    Pri_Nstrands = repmat((MinPriNstrands(1):1:MaxPriNstrands(1))', length(MaxPriNstrands), 1);

    % Calculate core loss (W)
    Iprms = Imax ./ sqrt(2);
    ue = ui ./ (1 + ui .* airgap ./ Le);
    Bm = u0.*Np.*Imax./Le.*ue;
    %(changed from /ue to *ue)

    % Window area
    Wa = H .* W;
    % Core volume (m3)
    Vcore = Ve;
    % Core weight (g)
    densityRow = CoreDensity(matno_record);
    assert(numel(densityRow)==numel(Vcore),'Core density mapping mismatch');
    Wcore = Vcore.*densityRow;
    % Check core loss
    Pcore = CoreLossMultiple .* Vcore .* K1 .* fs.^alpha .* Bm.^beta;

    %% Wire

    % Recalculate ds
    skindepth = 1 ./ sqrt(pi * fs * u0 / rou);
    % Overwritten variable
    ds = max(skindepth, MinLitzDia * ones(size(skindepth))); %#ok<NASGU> % take the skin depth litz
    ds = MinLitzDia * ones(size(skindepth)); % take the smallest litz
    
    % Primary wire diameter (m)
    Pri_WireSize = sqrt(Pri_Nstrands.*pi.*ds.^2./4./LitzFactor./pi).*2;
    
    % Primary wire diameter (m) including the insulation layer
    Pri_FullWireSize = Pri_WireSize+((Vpri./Np)./dielectricstrength_insulation).*2;
    % changed from Vin to Vpri./Np

    % Packing factors
    CopperPacking = (pi.*Pri_WireSize.^2.*Np./4)./(H.*W);
    OverallPacking = (pi.*Pri_FullWireSize.^2.*Np./4)./(H.*W);
    
    % Winding structures
    % Core insulation thickness needed
    CoreInsulationThickness = Vinsulation_max./dielectricstrength_insulation;
    
    % Turns per layer
    Pri_PerLayer = floor(Np./Mlp);
    
    % Total length of windings
    % For XcoreCoreShapeIndex == 1, EE cores
    TLp = Np.*2.*(PriW + PriH + 4*CoreInsulationThickness + ...
        2.*Mlp.*Pri_FullWireSize);

    % Recalculate XcoreCoreShapeIndex == 2, ER cores
    SelecIndex = find(XcoreCoreShapeIndex == 2);
    TLp(SelecIndex) = 2.*pi.*Np(SelecIndex).*(PriW(SelecIndex)./2 + ...
        CoreInsulationThickness(SelecIndex) + 0.5.*Mlp(SelecIndex).* ...
        Pri_FullWireSize(SelecIndex));
    
    % Copper Loss
    PriKlayer = sqrt(pi.*Pri_Nstrands).*ds./2./(Pri_WireSize);
    Pri_xp = ds./2./skindepth.*sqrt(pi.*PriKlayer);
    Pri_Rdc = rou.*TLp./(pi.*Pri_WireSize.^2./4);
    Pri_Fr = Pri_xp.*((sinh(2.*Pri_xp) + sin(2.*Pri_xp))./(cosh(2.*Pri_xp) ...
        - cos(2.*Pri_xp)) + 2.*(Mlp.^2.*Pri_Nstrands - 1)./3.*...
        (sinh(Pri_xp) - sin(Pri_xp))./(cosh(Pri_xp) + cos(Pri_xp)));
    Pri_Rac = Pri_Rdc.*Pri_Fr;
    Pcopper = Iprms.^2.*Pri_Rac;
    
    % Calculate temperature rise
    %Rth = 16.31.*1e-3.*(Ac.*Wa).^(-0.405);
    Rth = 16.31.*((Ac*1e4).*(Wa*1e4)).^(-0.405);
    %(GPT says this is a unit mismatch and to replace with the 1e-3 version,
    % since 16.31 is in cm, might be unit mismatch)

    Tafterloss = Rth.*(Pcopper + Pcore) + 25;
    
    % Calculate the weight
    WeightPri_copper = pi * Pri_WireSize.^2./4.*TLp.*CopperDensity;
    WeightPri_Insu = pi.*(Pri_FullWireSize.^2 - Pri_WireSize.^2)./4.*TLp.*WireInsulationDensity;
    
    % For XcoreCoreShapeIndex == 1, EE cores
    WeightCore_Insu = (2 .* H .* (PriW + 2 * PriH) + 4 .* W .* (PriW + 2 * PriH) + ...
        H .* (2.* PriW + 2 * PriH)) .* CoreInsulationThickness .* CoreInsulationDensity;
    
    % Recalculate XcoreCoreShapeIndex == 2, ER cores
    SelecIndex = find(XcoreCoreShapeIndex == 2);
    WeightCore_Insu(SelecIndex) = (sqrt(2) * pi * H(SelecIndex) .* PriW(SelecIndex) + ...
        sqrt(2) * pi * 2 * W(SelecIndex) .* PriW(SelecIndex) + ...
        H(SelecIndex)*pi.*PriW(SelecIndex)).*CoreInsulationThickness(SelecIndex).*CoreInsulationDensity;
    
    TotalWeight = Wcore + WeightPri_copper + WeightPri_Insu+WeightCore_Insu;
    %(added weightpri_insu)

    % Filter the designs
    B_index = find(Bm <= BSAT * BSAT_discount);
    P_loss_index = find((Pcopper + Pcore) <= Po * (1 - etaInductor));
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);
    TotalWeight_index = find(TotalWeight <= MaxWeight);
    
    OverallPackingmin_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactor);
    Mlp_index = find(Mlp .* Pri_FullWireSize <= W - 2 * CoreInsulationThickness);
    Pri_PerLayer_index = find(Pri_PerLayer.*Pri_FullWireSize<H-2*CoreInsulationThickness);
    
    idxList = {B_index, P_loss_index, Tafterloss_index, Tmin_index, TotalWeight_index, ...
               OverallPackingmin_index, OverallPackingmax_index, Mlp_index, Pri_PerLayer_index};
    
    Index_Meet_All = idxList{1};
    for k = 2:numel(idxList)
        Index_Meet_All = intersect(Index_Meet_All, idxList{k});
    end

    % Sort by total weight and keep only the lightest one
    [~, SortIndex] = sort(TotalWeight(Index_Meet_All));
    if ~isempty(SortIndex)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1));

        Volume=(Wcore(TotalWeightSortIndex))./CoreDensity+...
            (WeightPri_copper(TotalWeightSortIndex))./CopperDensity+...
            ((WeightPri_Insu(TotalWeightSortIndex)+WeightCore_Insu ...
            (TotalWeightSortIndex)))./CoreInsulationDensity;

        % Preallocate Design_inductor
        Design_inductor = zeros(numel(TotalWeightSortIndex), 43);

        Design_inductor(:, 1)  = Po(TotalWeightSortIndex);
        Design_inductor(:, 2)  = Vin(TotalWeightSortIndex);
        Design_inductor(:, 3)  = Vpri(TotalWeightSortIndex);
        Design_inductor(:, 4)  = Vinsulation_max(TotalWeightSortIndex);
        Design_inductor(:, 5)  = fs(TotalWeightSortIndex);
        Design_inductor(:, 6)  = matno_record(TotalWeightSortIndex);
        Design_inductor(:, 7)  = matfs(TotalWeightSortIndex);
        Design_inductor(:, 8)  = PriW(TotalWeightSortIndex);
        Design_inductor(:, 9)  = PriH(TotalWeightSortIndex);
        Design_inductor(:,10)  = Ac(TotalWeightSortIndex);
        Design_inductor(:,11)  = H(TotalWeightSortIndex);
        Design_inductor(:,12)  = W(TotalWeightSortIndex);
        Design_inductor(:,13)  = Np(TotalWeightSortIndex);
        Design_inductor(:,14)  = Bm(TotalWeightSortIndex);
        Design_inductor(:,15)  = Pri_WireSize(TotalWeightSortIndex);
        Design_inductor(:,16)  = Pri_FullWireSize(TotalWeightSortIndex);
        Design_inductor(:,17) = Imax(TotalWeightSortIndex) ./ ...
            (pi * Pri_Nstrands(TotalWeightSortIndex) .* ds(TotalWeightSortIndex).^2 / 4);
        Design_inductor(:,18)  = Pri_Nstrands(TotalWeightSortIndex);
        Design_inductor(:,19)  = Pri_PerLayer(TotalWeightSortIndex);
        Design_inductor(:,20)  = Mlp(TotalWeightSortIndex);
        Design_inductor(:,21)  = CopperPacking(TotalWeightSortIndex);
        Design_inductor(:,22)  = OverallPacking(TotalWeightSortIndex);
        Design_inductor(:,23)  = Pcore(TotalWeightSortIndex);
        Design_inductor(:,24)  = Pcopper(TotalWeightSortIndex);
        Design_inductor(:,25)  = Wcore(TotalWeightSortIndex);
        Design_inductor(:,26)  = WeightPri_copper(TotalWeightSortIndex);
        Design_inductor(:,27)  = WeightPri_Insu(TotalWeightSortIndex);
        Design_inductor(:,28)  = WeightCore_Insu(TotalWeightSortIndex);
        Design_inductor(:,29)  = TotalWeight(TotalWeightSortIndex);
        Design_inductor(:,30)  = Tafterloss(TotalWeightSortIndex);
        Design_inductor(:,31)  = L(TotalWeightSortIndex);
        Design_inductor(:,32)  = airgap(TotalWeightSortIndex);
        Design_inductor(:,33)  = XcoreIndex(TotalWeightSortIndex);
        Design_inductor(:,34)  = LCC_Q  * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,35)  = LCC_f0 * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,36)  = LCC_A  * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,37)  = LCC_K  * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,38)  = LCC_RT * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,39)  = LCC_Ls * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,40)  = LCC_Cs * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,41)  = LCC_Cp * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,42)  = LCC_GT * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,43) = Volume;
        y = Design_inductor;
    else
        y = zeros(1,43);
        warning('No successful results Test2, Inductor');
    end
end
%toc
end

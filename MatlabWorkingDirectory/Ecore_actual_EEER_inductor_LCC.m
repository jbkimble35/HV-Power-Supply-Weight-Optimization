function y = Ecore_actual_EEER_inductor_LCC(raw,raw1,raw2,raw3,raw4,raw5,raw6, ...
    Vin_range, G_range, Po_range, fs_range, Ls_range, Imax_range, Winding_Pattern, ...
    LCC_Q, LCC_f0, LCC_A, LCC_K, LCC_RT, LCC_Ls, LCC_Cs, LCC_Cp, LCC_GT)

% Tunable Parameters
%% -------------------------------------------------------------------------------------

% Insulation
%-------------------------------------------

% Need to add interlayer tape
layerTapeUse = true;
enamelThickness = 60e-6;
kaptonDielStrength = 0.5*200e5;
kaptonThickness = 60e-6;
MinTapeMargin = 5e-4;
kaptonDensity = 1.42e6; % g/m^3
cover_fac  = 1.05; % overlap of tape on other areas

% Dielectric strength of the insulation material (V/m), discount 50%
dielectricstrength_insulation = 0.5 * 200e5; % TEFLON
% g/m^3, density of core insulation materials
CoreInsulationDensity = 2.2e6; % TEFLON
% g/m^3, density of wire insulation materials
WireInsulationDensity = 2.2e6; % TEFLON

% Inductor parameters
%-------------------------------------------

% Lowest allowed inductor efficiency
etaInductor = 0.90;
% Max allowable temperature (C)
Tmax = 100;
% Min allowable temperature (C)
Tmin = 25;
% Maximum allowable weight (g)
MaxWeight = 10000;
% Air gap (m)
mingap = 0;

% Winding and Wire Parameters
%------------------------------------------

% Minimum turns
MinWinding = 1;
% Maximum turns
MaxWinding = 200;
% Incremental winding
IncreN = 1;
% Maximum layer of winding
MaxMl = 10;
% Incremental layers. The layers of a transformer reference each wrap of
% turns that fills the window height before moving on to the next level.
% Once one layer fills, the next layer is wound on top, seperated by an
% insulation layer.
IncreMl = 1;
% Minimal wire diameter (m)
MinWireDia = 0.079/1000; % AWG28, 0.35 mm is AWG29, 0.079 is AWG40
% Max allowable current density in the wire (A/m^2)
% 500A/cm^2 is the upper bound recommended, but without active cooling, and
% since the magnetics are thermally insulated, less is assumed
Jwmax = 3e6;
% Minimal litz diameter one can get (m)
MinLitzDia = 0.05024 / 1000; % AWG44 % 0.0316 is AWG48, 0.03983 is AWG46
% g/m^3, density of copper
CopperDensity = 8.96*1000*1000;
% Copper wire multiple to reduce resistive losses
CuMult = 1;

% Discount factors
%----------------------------------------

% Bmax discount factor
BSAT_discount = 0.75;
% Actual core loss is always higher than the calculated
CoreLossMultiple = 1;
% Maximum packing factor (copper area compared with total window area)
maxpackingfactor = 0.7;
% Minimum packing factor
minpackingfactor = 0.01;
% Winding factor of litz wire, assuming only 80% of wire size is copper
LitzFactor = 0.8;

% Electrical Parameters
%----------------------------------------

% Electrical constants. Normally there is no need to change
% ohm*m, resistivity of copper at 100C
rou = 2.3*1e-8;
% H/A·m^2, permeability of free space
u0 = 4*pi*10^(-7);


% Function Body
%% -------------------------------------------------------------------------------------

% Parse core loss material maps from CoreLossData.xlsx
[m1,n1]   = size(raw1);
LCoreFreq = cell2mat(raw1(2:m1,3:n1));
[m1,n1]   = size(raw2);
LCoreBfield = cell2mat(raw2(2:m1,3:n1));
[m1,n1]   = size(raw3);
LCorePloss   = cell2mat(raw3(2:m1,3:n1));
% Only parses column C
[m1,~]    = size(raw4);
LCoreBSAT = cell2mat(raw4(2:m1,3));
% Only parses column C
[m1,~]    = size(raw5);
LCoreInitialMU   = cell2mat(raw5(2:m1,3));

[m1,~]    = size(raw6);
CoreDensity   = cell2mat(raw6(2:m1,3))*1000000;

% Reference loss level and PF-factor for Steinmetz equations
Pbar = 500;      % mW/cm^3
PFfactor = 1;

% Build Steinmetz parameter sets around the target frequency
% ----------------------------------------------

% Number of materials excluding the header column
NoMat = m1-1;
% Preallocates flag that will be 0 if the material doesn't have a frequency
% near the design point, and 1 if it does, for each of the material indexes
FreqFlag = zeros(size(1:1:NoMat));
% Sweeps all of the materials. This is the core loss curve fitting loop for
% each material. It processes the datasheet parameters into Steinmetz parameters.
for i = 1:1:NoMat
    % Obtains all frequency values for the material, ignoring NaN values.
    DataSheetFreq = LCoreFreq(i,~isnan(LCoreFreq(i,:)));
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
        ConstantA(i,j) = (log10(LCorePloss(i,2*j)) - log10(LCorePloss(i,2*j-1)))/( ...
            log10(LCoreBfield(i,2*j)) - log10(LCoreBfield(i,2*j-1)));
        ConstantB(i,j) = log10(LCorePloss(i,2*j)) - ConstantA(i ,j)*log10(LCoreBfield( ...
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
            beta_range(i,j) = log10(LCorePloss(i,2*j)/LCorePloss(i,2*j-1))/log10( ...
                LCoreBfield(i,2*j)/LCoreBfield(i,2*j-1));
            % Extrapolate third point on previous line (same B2)
            XCorePloss_3rd(i,j) = 10.^(ConstantA(i,j-1)*log10(LCoreBfield(i,2*j))+ ...
                ConstantB(i,j-1));
            % alpha exponent from frequency dependence
            alpha_range(i,j) = log10(XCorePloss_3rd(i,j)/LCorePloss(i,2*j))/log10( ...
                DataSheetFreq(2*j-3)/DataSheetFreq(2*j-1)); %(f2/f1)^alpha = P2/P1;
            % K coefficient in Steinmetz equation
            K1_range(i,j) = LCorePloss(i,2*j)/(LCoreBfield(i,2*j)^beta_range(i,j))/( ...
                DataSheetFreq(2*j-1)^alpha_range(i,j)); %mW/cm3

            % Populate j-1 if missing
            if j == 2
                beta_range(i,j-1)  = beta_range(i,j);
                alpha_range(i,j-1) = alpha_range(i,j);
                K1_range(i,j-1) = LCorePloss(i,2*j-2)/(LCoreBfield(i,2*j-2)^ ...
                    beta_range(i,j-1))/(DataSheetFreq(2*j-3)^alpha_range(i,j-1));            
            end
        end
    end
end


% Core size
%-------------------------------------------

[m1,~] = size(raw);
InductorCoreIndex = cell2mat(raw(2:m1,1));
% Effective core volume in m^3
LcoreVe = cell2mat(raw(2:m1,3))/(1000^3);
% Cross-sectional area in m^2
LcoreAe = cell2mat(raw(2:m1,4))/(1000^2);
% Magnetic path length in m
LcoreLe = cell2mat(raw(2:m1,5))/1000;

% 1 for EE, 2 for ER, will be 3 for U and 4 for UR
LcoreCoreShapeIndex = cell2mat(raw(2:m1,6));
% Primary winding window dimensions in m
LcorePriW = cell2mat(raw(2:m1,8))/1000;
LcorePriH = cell2mat(raw(2:m1,9))/1000;

LcoreWindowW = cell2mat(raw(2:m1,12))/1000;
LcoreWindowH = cell2mat(raw(2:m1,13))/1000;
ShuffleIndex = 1:1:length(InductorCoreIndex);

% DESIGN SWEEP
%% ------------------------------------------------------------------------

% Material indices with data near fs
CoreMatIndexSweep = find(FreqFlag);  

% Need to make this more efficient somehow. Currently eats up too much RAM
% Creates 10-D array of all ranges of values for sweeping.

[Po, fs, Vin, G, Ls, Imax, matno_record, ShuffleLcoreIndex, Np, Mlp] = ndgrid( ...
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
ShuffleLcoreIndex = reshape(ShuffleLcoreIndex,[],1);
Np = reshape(Np,[],1);
Mlp = reshape(Mlp,[],1);

% Maximum voltage stress for each candidate
Vinsulation_max = Vin.*G;

% Finds BSAT and relative permittivity for each candidate design, using matno_record
% as the index deciding which material is assigned to which index.
ui = LCoreInitialMU(matno_record);

BSAT = LCoreBSAT(matno_record);

% ShuffleLcoreIndex is a very large column vector containing values of indices
% for rows in the CoreSizeData file, and in the following lines, creates very large column vectors
% corresponding to the dimensions for each core geometry chosen. Because
% each core geometry consists of 8 parameters, this is just unfolding those.
Ve      = LcoreVe(ShuffleLcoreIndex);
W       = LcoreWindowW(ShuffleLcoreIndex);
H       = LcoreWindowH(ShuffleLcoreIndex);
PriW    = LcorePriW(ShuffleLcoreIndex);
PriH    = LcorePriH(ShuffleLcoreIndex);
Ac      = LcoreAe(ShuffleLcoreIndex);
Le      = LcoreLe(ShuffleLcoreIndex);
LcoreIndex = InductorCoreIndex(ShuffleLcoreIndex);
LcoreCoreShapeIndex = LcoreCoreShapeIndex(ShuffleLcoreIndex);


% Primary winding voltage seen by inductor
Vpri    = Vin.*G;
% Shorthand for ease of readability
L       = Ls;
% Effective air gap length based on the magnetic circuit relation between
% inductance, turns, effective length, cross-sectional area, etc.

airgap = u0.*Ac.*Np.^2./L-Le./ui;
% Eliminate unrealistic elements: require a feasible min gap
KeepAirGap = intersect(find(airgap>=mingap), find(airgap<=0.2*Le));

% Safety check for air gaps. If the computed air gap is negative, it means
% the inductance isn't enough.
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
%L_calc = u0.*Np.^2.*Ac.*ui./(Le + ui.*airgap);
%L_index = find(abs(L_calc - L)./L <= 0.15);
%KeepIndex = intersect(KeepIndex, L_index);

% Debug
if isempty(KeepIndex)
    warning('No Successful Indexes');
    y = zeros(1,37);
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
LcoreIndex         = LcoreIndex(KeepIndex);
LcoreCoreShapeIndex= LcoreCoreShapeIndex(KeepIndex);

% Filters out <=0 core loss at all frequencies, material frequencies >40% away from fs,
% and updates candidate frequency points into matfs.
FsnoNonzero = F_atPv_500(matno_record,:)>0;
FsnoIndex = abs(fs-F_atPv_500(matno_record,:))./fs<=0.4;
matfsIndex = FsnoNonzero.*FsnoIndex;
matfs = F_atPv_500(matno_record,:).*matfsIndex;

% Converts K1 from mW/cm3 to W/m3 and builds core-loss steinmetz parameters
% K1, alpha, and beta for the frequencies given by matfs. Then rowIdcs
% stores the indices in the arrays that correspond to valid frequency
% points.
K1 = K1_range(matno_record,:).*matfsIndex*1000;
alpha = alpha_range(matno_record,:).*matfsIndex;
beta = beta_range(matno_record,:).*matfsIndex;
[rowIdcs, ~] = find(matfs>0);

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
LcoreIndex          = repelem(LcoreIndex(UniqueRowIdcs), ColDuplicate);
LcoreCoreShapeIndex = repelem(LcoreCoreShapeIndex(UniqueRowIdcs), ColDuplicate);

% Reformat loss data into one non-zero vector
matfs = nonzeros(reshape(matfs(UniqueRowIdcs,:)',[],1));
K1 = nonzeros(reshape(K1(UniqueRowIdcs,:)',[],1));
beta = nonzeros(reshape(beta(UniqueRowIdcs,:)',[],1));
alpha = nonzeros(reshape(alpha(UniqueRowIdcs,:)',[],1));


% For the remaining indexes, compute the more detailed parameters 
% like losses, size, subsequent weight
% -------------------------------------------------------------------

if (isempty(Po))
    warning('No successful inductor results Test1, Inductor.');
    y = 0;
else


    % Flux and core loss
    % -------------------------------------

    ue = ui./(1+ui.*airgap./Le);
    Bm = u0.*Np.*Imax./Le.*ue;

    % Core weight (g)
    Wcore = Ve.*CoreDensity(matno_record);

    % Check core loss
    Pcore = CoreLossMultiple.*Ve.*K1.*fs.^alpha.*Bm.^beta;

    % Preallocate
    TLp = zeros(size(Np));

    % Determine wire type, size, and num of strands if Litz
    % -----------------------------------------------

    % RMS current
    Iprms=Imax./sqrt(2);
    % AC skin depth
    skindepth=1./sqrt(pi.*fs.*u0./rou);
    % Area required of wire m^2
    Areq_p=CuMult.*Iprms./Jwmax;
    % solid equivalent diameter
    dsolid=2.*sqrt(Areq_p./pi);
    % Solid vs. litz
    useSolid=dsolid<=2.*skindepth;
    % Litz diameter
    dstrand_litz=max(MinLitzDia,2.*skindepth);
    % strand cross section area
    Astrand=pi.*(dstrand_litz./2).^2;
    % Number of strands default to 1
    Pri_Nstrands=ones(size(Iprms));
    % Where not solid, use litz number of strands
    Pri_Nstrands(~useSolid)=ceil(Areq_p(~useSolid)./Astrand(~useSolid));
    % use solid diameter if solid, bundle diameter if litz
    Pri_WireDia=max(MinWireDia,dsolid);
    idLitz=~useSolid;
    if any(idLitz)
        Pri_WireDia(idLitz)=2.*sqrt(Pri_Nstrands(idLitz).*Astrand(idLitz).*LitzFactor./pi);
    end
    
    % Strand diameter
    Pri_ds=max(MinWireDia,dsolid);
    Pri_ds(idLitz)=dstrand_litz(idLitz);
    % Full wire size with insulation
    Pri_FullWireDia=Pri_WireDia+2.*Vpri./dielectricstrength_insulation;
    if layerTapeUse
        Pri_FullWireDia = Pri_WireDia + enamelThickness.*2;
    end

    A_pri_cu=(pi.*(Pri_WireDia.^2))./4;
    A_pri_full=(pi.*(Pri_FullWireDia.^2))./4;

    CopperPacking=(A_pri_cu.*Np)./(H.*W);
    OverallPacking=(A_pri_full.*Np)./(H.*W);

    CoreInsulationThickness=Vpri./dielectricstrength_insulation;
    
    % Computes mean length of turn for pri and sec, accounting for geometry
    % and winding pattern
    %-----------------------------------------------------------------------------

    % Turns per layer
    Pri_PerLayer = floor(Np./Mlp);

    if layerTapeUse
        numTapePerLayerPri = ceil((Vpri./Mlp)./(kaptonDielStrength*kaptonThickness));
        tTapePri = max(Mlp-1,0).*numTapePerLayerPri.*kaptonThickness;
    else
        numTapePerLayerPri = zeros(size(Mlp));
        tTapePri = zeros(size(Mlp));
    end

    isEE=(LcoreCoreShapeIndex==1);
    isER=(LcoreCoreShapeIndex==2);
    isU=(LcoreCoreShapeIndex==3);
    isUR=(LcoreCoreShapeIndex==4);
    
    if Winding_Pattern==1
        TLp(isEE|isU)=2.*Np(isEE|isU).*(PriW(isEE|isU)+PriH(isEE|isU)+ ...
            4.*CoreInsulationThickness(isEE|isU)+2.*tTapePri(isEE|isU)+ ...
            2.*Mlp(isEE|isU).*Pri_FullWireDia(isEE|isU));
        TLp(isER|isUR)=2.*pi.*Np(isER|isUR).*(PriW(isER|isUR)./2+ ...
            CoreInsulationThickness(isER|isUR)+0.5.*tTapePri(isER|isUR)+ ...
            0.5.*Mlp(isER|isUR).*Pri_FullWireDia(isER|isUR));
    end

    if Winding_Pattern==2
        TLp(isEE|isU)=2.*Np(isEE|isU).*(PriW(isEE|isU)+PriH(isEE|isU)+ ...
            4.*CoreInsulationThickness(isEE|isU)+2.*tTapePri(isEE|isU)+ ...
            2.*Mlp(isEE|isU).*Pri_FullWireDia(isEE|isU));
        TLp(isER|isUR)=2.*pi.*Np(isER|isUR).*(PriW(isER|isUR)./2+ ...
            CoreInsulationThickness(isER|isUR)+0.5.*tTapePri(isER|isUR)+ ...
            0.5.*Mlp(isER|isUR).*Pri_FullWireDia(isER|isUR));
    end

    Mlp_index=find((Mlp.*Pri_FullWireDia)<=(W-2.*CoreInsulationThickness));
    Pri_PerLayer_index=find((Pri_PerLayer.*Pri_FullWireDia)<=(H-2.*CoreInsulationThickness));
    
    % Calculate Copper Loss & Temp rise
    %--------------------------------------------------------

    PriKlayer=sqrt(pi.*Pri_Nstrands).*Pri_ds./(2.*Pri_WireDia);
    Pri_xp=Pri_ds./(sqrt(pi.*PriKlayer).*2.*skindepth);

    % Fixed overestimation by litzfactor
    Pri_Rdc = rou.*TLp./(pi.*dsolid.^2./4);
    Pri_Fr=Pri_xp.*((sinh(2.*Pri_xp)+sin(2.*Pri_xp))./(cosh(2.*Pri_xp)-cos(2.*Pri_xp))+2.* ...
        (Mlp.^2.*Pri_Nstrands-1)./3.*(sinh(Pri_xp)-sin(Pri_xp))./(cosh(Pri_xp)+cos(Pri_xp)));
    Pri_Rac=Pri_Rdc.*Pri_Fr;
    Pcopper=Iprms.^2.*Pri_Rac;

    Wa=2.*H.*W.*(isEE|isER)+H.*W.*(isU|isUR);
    Rth=16.31e-3.*(Ac.*Wa).^(-0.405);
    Tafterloss=Rth.*(Pcopper+Pcore)+25;

    % Calculate the weight of core insulation
    %----------------------------------------------------------------
    
    % Must add weight of potting and weight of tape here, this is just
    % core liner of teflon

    WeightCore_Insu=zeros(size(Np));
    
    % EE
    WeightCore_Insu(isEE)=( ...
        2.*H(isEE).*(PriW(isEE)+2.*PriH(isEE)) + ...
        4.*W(isEE).*(PriW(isEE)+2.*PriH(isEE)) + ...
        H(isEE).*(2.*PriW(isEE)+2.*PriH(isEE)) ) ...
        .*CoreInsulationThickness(isEE).*CoreInsulationDensity;
    
    % ER
    WeightCore_Insu(isER)=( ...
        sqrt(2).*pi.*H(isER).*PriW(isER) + ...
        sqrt(2).*pi.*2.*W(isER).*PriW(isER) + ...
        H(isER).*pi.*PriW(isER) ) ...
        .*CoreInsulationThickness(isER).*CoreInsulationDensity;
    
    % U
    WeightCore_Insu(isU)=( ...
        2.*H(isU).*(PriW(isU)+2.*PriH(isU)) + ...
        2.*W(isU).*(PriW(isU)+2.*PriH(isU)) + ...
        H(isU).*(2.*PriW(isU)+2.*PriH(isU)) ) ...
        .*CoreInsulationThickness(isU).*CoreInsulationDensity;
    
    % UR
    WeightCore_Insu(isUR)=( ...
        2.*H(isUR).*(PriW(isUR)+2.*PriH(isUR)) + ...
        pi.*W(isUR).*PriW(isUR) + ...
        H(isUR).*pi.*PriW(isUR) ) ...
        .*CoreInsulationThickness(isUR).*CoreInsulationDensity;

    % Copper & wire-insulation weights
    %---------------------------------------------

    WeightPri_copper = pi.*Pri_WireDia.^2./4.*TLp.*CopperDensity;
    WeightPri_Insu=(pi.*(Pri_FullWireDia.^2-Pri_WireDia.^2)./4).*TLp.*WireInsulationDensity;

    % Interlayer Tape Weight
    %---------------------------------------------

    if ~layerTapeUse
        Weight_InterlayerTape = zeros(size(Mlp));  % g
        V_tape = zeros(size(Mlp)); % g
        tapeMargin = 0;
    else
        tapeMargin = max(0.02*H, MinTapeMargin);
        % total build incl. tape (for average circumference)
        a1 = Mlp.*Pri_FullWireDia + max(Mlp-1,0).*numTapePerLayerPri.*kaptonThickness; % m
    
        % base per-turn length near inner radius
        Lbase_p(isEE|isU)  = 2.*(PriW(isEE|isU) + PriH(isEE|isU) + 4.*CoreInsulationThickness(isEE|isU));
        Lbase_p(isER|isUR) = 2.*pi.*(PriW(isER|isUR)./2 + CoreInsulationThickness(isER|isUR));
        Lbase_p = Lbase_p';

        % average circumference of tape wraps
        Lavg_il_p = Lbase_p + pi.*(a1./2);
    
        % number of inter-layer boundaries
        nBndPri = max(Mlp-1,0);
    
        % total tape length (all boundaries × plies × overlap factor)
        L_tape_total = cover_fac.*(nBndPri.*numTapePerLayerPri.*Lavg_il_p); % m
    
        % tape width and volume
        w_tape = H+2.*tapeMargin;                 % m
        V_tape = kaptonThickness.*w_tape.*L_tape_total;  % m^3
    
        % mass in grams (KaptonDensity in g/m^3)
        Weight_InterlayerTape = kaptonDensity.*V_tape;  % g
    end
        
    % Total weight
    %---------------------------------------------

    TotalWeight=Wcore+WeightPri_copper+WeightPri_Insu+WeightCore_Insu+Weight_InterlayerTape;

    % Filter the designs
    %---------------------------------------------

    % Filter by B
    B_index = find(Bm < BSAT*BSAT_discount);
    [Bmin,BminIndex] = min(Bm);
        
    % Filter by Temperature and Power Loss
    P_loss_index = find(Pcopper + Pcore <= Po./etaInductor - Po);
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);
    [Tminimum,TminValIndex] = min(Tafterloss);
    [Pmin,PminValIndex] = min(Pcopper+Pcore);

    % Filter by weight
    TotalWeight_index = find(TotalWeight <= MaxWeight);
    [WMin,WminValIndex] = min(TotalWeight);

    % Filter by packing factor min and max
    if layerTapeUse
        % might be incorrect due to approx. as H*t 
        if Winding_Pattern == 2
            CopperPacking = ((pi.*(Pri_WireDia.^2)./4).*Np) ./ ((H-2.*tapeMargin).*W);
            OverallPacking = ( ((pi.*(Pri_FullWireDia.^2)./4).*Np) ...
                + (H .* ( max(Mlp-1,0).*numTapePerLayerPri.*kaptonThickness )) ) ...
                ./ ((H-2.*tapeMargin).*W);
        else 
            CopperPacking = ((pi.*(Pri_WireDia.^2)./4).*Np) ./ ((H-2.*tapeMargin).*W);
            OverallPacking = ((pi.*(Pri_FullWireDia.^2)./4).*Np) ./ ((H-2.*tapeMargin).*W);
        end
    end
    OverallPackingmin_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactor);
    [PackingMin,PackingMinValIndex] = max(OverallPacking);
    [PackingMax,PackingMaxValIndex] = min(OverallPacking);

    % Filter by layer width and height
    [MlpMax,MlpMaxValIndex] = max(Mlp.*Pri_FullWireDia);
    [PriPerLayerMax,PPLMaxIndex] = max(Pri_PerLayer.*Pri_FullWireDia);

    if isempty(P_loss_index)
        fprintf("Inductor Bottlenecked. Min Power loss out of all candidates: %.2f Index: %d",Pmin,PminValIndex);
        y = zeros(1,37);
        return
    end
    if isempty(Tafterloss_index)
        fprintf("Inductor Bottlenecked. Min T out of all candidates: %.2f Index: %d",Tminimum,TminValIndex);
        y = zeros(1,37);
        return
    end
    if isempty(B_index)
        fprintf("Inductor Bottlenecked. Min B out of all candidates: %.2f Index: %d",Bmin,BminIndex);
        y = zeros(1,37);
        return
    end
    if isempty(TotalWeight_index)
        fprintf("Inductor Bottlenecked. Min Weight out of all candidates: %.2f Index: %d",WMin,WminValIndex);
        y = zeros(1,37);
        return
    end
    if isempty(OverallPackingmin_index)
        fprintf("Inductor Bottlenecked. Min packing factor out of all candidates: %.2f Index: %d",PackingMin,PackingMinValIndex);
        y = zeros(1,37);
        return
    end
    if isempty(OverallPackingmax_index)
        fprintf("Inductor Bottlenecked. Max packing factor out of all candidates: %.2f Index: %d",PackingMax,PackingMaxValIndex);
        y = zeros(1,37);
        return
    end
    if isempty(Mlp_index)
        fprintf("Inductor Bottlenecked. Max layers width of all candidates: %.2f Index: %d",MlpMax,MlpMaxValIndex);
        y = zeros(1,37);
        return
    end
    if isempty(Pri_PerLayer_index)
        fprintf("Inductor Bottlenecked. Max layer height out of all candidates: %.2f Index: %d",PriPerLayerMax,PPLMaxIndex);
        y = zeros(1,37);
        return
    end

    Index_Meet_All = B_index;
    Index_Meet_All = intersect(Index_Meet_All, P_loss_index);
    Index_Meet_All = intersect(Index_Meet_All, Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All, Tmin_index);
    Index_Meet_All = intersect(Index_Meet_All, TotalWeight_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmin_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmax_index);
    Index_Meet_All = intersect(Index_Meet_All, Mlp_index);
    Index_Meet_All = intersect(Index_Meet_All, Pri_PerLayer_index);


    % Build Results Table 
    %% -----------------------------------------------------------------------------

    % Keep only the lightest design
    [~, SortIndex] = sort(TotalWeight(Index_Meet_All));
    if ~isempty(SortIndex)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1));
        Volume = Ve(TotalWeightSortIndex) ...
           + WeightPri_copper(TotalWeightSortIndex)./CopperDensity ...
           + WeightPri_Insu(TotalWeightSortIndex)./WireInsulationDensity ...
           + WeightCore_Insu(TotalWeightSortIndex)./CoreInsulationDensity;

        % Preallocate Design_inductor
        Design_inductor = zeros(numel(TotalWeightSortIndex), 37);

        Design_inductor(:, 1)  = Po(TotalWeightSortIndex);
        Design_inductor(:, 2)  = Vin(TotalWeightSortIndex);
        Design_inductor(:, 3)  = Vpri(TotalWeightSortIndex);
        Design_inductor(:, 4)  = fs(TotalWeightSortIndex);
        Design_inductor(:, 5)  = matno_record(TotalWeightSortIndex);
        Design_inductor(:, 6)  = matfs(TotalWeightSortIndex);
        Design_inductor(:, 7)  = Np(TotalWeightSortIndex);
        Design_inductor(:, 8)  = Bm(TotalWeightSortIndex);
        Design_inductor(:, 9)  = Pri_WireDia(TotalWeightSortIndex);
        Design_inductor(:,10)  = Pri_FullWireDia(TotalWeightSortIndex);
        Design_inductor(:,11)  = Imax(TotalWeightSortIndex) ./ ...
            (pi * Pri_Nstrands(TotalWeightSortIndex) .* Pri_ds(TotalWeightSortIndex).^2 / 4);
        Design_inductor(:,12)  = Pri_Nstrands(TotalWeightSortIndex);
        Design_inductor(:,13)  = Pri_PerLayer(TotalWeightSortIndex);
        Design_inductor(:,14)  = Mlp(TotalWeightSortIndex);
        Design_inductor(:,15)  = CopperPacking(TotalWeightSortIndex);
        Design_inductor(:,16)  = OverallPacking(TotalWeightSortIndex);
        Design_inductor(:,17)  = Pcore(TotalWeightSortIndex);
        Design_inductor(:,18)  = Pcopper(TotalWeightSortIndex);
        Design_inductor(:,19)  = Wcore(TotalWeightSortIndex);
        Design_inductor(:,20)  = WeightPri_copper(TotalWeightSortIndex);
        Design_inductor(:,21)  = WeightPri_Insu(TotalWeightSortIndex);
        Design_inductor(:,22)  = WeightCore_Insu(TotalWeightSortIndex);
        Design_inductor(:,23)  = TotalWeight(TotalWeightSortIndex);
        Design_inductor(:,24)  = Tafterloss(TotalWeightSortIndex);
        Design_inductor(:,25)  = L(TotalWeightSortIndex);
        Design_inductor(:,26)  = airgap(TotalWeightSortIndex);
        Design_inductor(:,27)  = LcoreIndex(TotalWeightSortIndex);
        Design_inductor(:,28)  = LCC_Q  * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,29)  = LCC_f0 * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,30)  = LCC_A  * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,31)  = LCC_K  * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,32)  = LCC_RT * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,33)  = LCC_Ls * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,34)  = LCC_Cs * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,35)  = LCC_Cp * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,36)  = LCC_GT * ones(numel(TotalWeightSortIndex),1);
        Design_inductor(:,37)  = Volume;

        y = Design_inductor;
    else
        y = zeros(1,37);
        warning('No successful results Test2, Inductor');
    end
end
%toc
end

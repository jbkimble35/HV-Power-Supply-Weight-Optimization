function y = Ecore_actual_EEER_xfmer_LCC(raw,raw1,raw2,raw3,raw4,raw5,raw6, ...
    Vppeak_range, Vspeak_range, Po_range, fs_range, Vinsulation_max_range, Winding_Pattern)

% Tunable Parameters
%% -------------------------------------------------------------------------------------

% The inductor is much lower voltage since it's on the primary side, so it
% doesn't need some of these insulation parameters.

% Insulation
%-------------------------------------------

% Potting spec.s
UsePotting      = true;
Potting_DS      = 20e6;      % V/m
CoronaMargin    = 2.0;       % Breakdown voltage safety factor

% Enamel/Insulation one-sided thickness
t_turn_p = 40e-6; %m 
t_turn_s = 60e-6; %m

% Interlayer tape (e.g., Kapton)
Tape_Interlayer_Thickness = 50e-6;    % m, per ply
Tape_Interlayer_DS        = 200e6;    % V/m, (≈200 kV/mm for Kapton)
Tape_Interlayer_Wraps     = 2;

% Thickness of epoxy between layers
Interlayer_Epoxy_Thickness = 0;       % m

% GPU computing options
%-------------------------------------------

useGPU = true;
useSingleOnGPU = false;

% Transformer parameters
%-------------------------------------------

% Minimum transformer efficiency
etaXfmer = 0.85;
% Max operating temp in Celsius
Tmax = 100;
% Min operating temp in Celsius
Tmin = 25;
% Max allowable current density in wire, (A/m^2)
Jwmax = 500*100*100;
% Minimal litz diameter (A/m^2)
MinLitzDia = 0.07874e-3; % AWG40
% Dielectric strength of insulation material (V/m) 25% derated
dielectricstrength_insulation = 0.75 * 200e3 * 100;
% Minimum primary windings
MinPriWinding = 1;
% Maximum primary windings
MaxPriWinding = 501;
% Incremental primary winding
IncreNp = 10;
% Maximum layer of primary winding
MaxMlp = 10;
% Incremental layer of primary winding
IncreMlp = 1;
% Maximum layer of secondary winding
MaxMls = 50;
% Incremental layer of secondary winding
IncreMls = 2;
% Minimum secondary wire diameter (m)
MinSecWireSize = 0.4e-3;             %#ok<NASGU>
% Max allowable transformer weight (g)
MaxWeight = 300000;
CoreInsulationDensity = 2.2e6;       % g/m^3 (Teflon)
WireInsulationDensity = 2.2e6;       % g/m^3 (Teflon)
% Saturation flux density derating
BSAT_discount = 0.9;
% Core loss multiplier
CoreLossMultiple = 1.5;
maxpackingfactor = 0.99;
minpackingfactor = 0.01;
% Litz copper fill
LitzFactor = 0.8;
% For solid core:
SolidGate = 1.5; % factor of skindepth threshold to go from solid to litz

% Electrical constants
%-------------------------------------------

% density of copper (g/m^3)
CopperDensity = 8.96e6;
% ohm*m, resistivity of copper at 100C
rou = 2.3e-8;
% conductivity of copper, inverse of rou
sigma = 1/rou; %#ok<NASGU>
% permittivity of free space HA/m^2
u0 = 4*pi*1e-7;
% permittivity of free space F/m
ebs10 = 8.854*1e-12; %#ok<NASGU>


% Body of function
%% --------------------------------------------------------------------------------------


if useSingleOnGPU
    toGPU    = @(x) gpuArray(single(x));
    toScalar = @(x) single(x);
else
    toGPU    = @(x) gpuArray(x);
    toScalar = @(x) x;
end
toHost = @(x) gather(x);

% Preallocate
Design = zeros(1,44);

[m1,n1] = size(raw1);
XCoreMAT = raw1(2:m1,2); %#ok<NASGU>
XCoreFreq = cell2mat(raw1(2:m1,3:n1));
[m1,n1] = size(raw2);
XCoreBfield    = cell2mat(raw2(2:m1,3:n1));
[m1,n1] = size(raw3);
XCorePloss   = cell2mat(raw3(2:m1,3:n1));

[NoMat,~] = size(XCoreFreq);
for i = 1:NoMat
    frow = XCoreFreq(i,~isnan(XCoreFreq(i,:)));
    brow = XCoreBfield(i,   ~isnan(XCoreBfield(i,:)));
    prow = XCorePloss(i,  ~isnan(XCorePloss(i,:)));

    if mod(numel(frow),2) ~= 0
        error('Data:OddFreqPairs','Material row %d has an odd count of frequency entries; must be pairs.', i);
    end
    if numel(brow) ~= numel(prow)
        error('Data:BPMismatch','Material row %d has %d B entries but %d Ploss entries.', i, numel(brow), numel(prow));
    end
    if numel(brow) ~= numel(frow)
        error('Data:FreqBPLength','Material row %d: Freq cols (%d) must equal B/Ploss cols (%d).', ...
              i, numel(frow), numel(brow));
    end
end

% Made BSAT, Density and MU read only 1 value per material.
[m1,~] = size(raw4);
XCoreBSAT  = cell2mat(raw4(2:m1,3));
[m1,~] = size(raw5);
XCoreMU    = cell2mat(raw5(2:m1,3));
[m1,~]    = size(raw6);
CoreDensity   = cell2mat(raw6(2:m1,3))*1000000;

Pbar = 500;      % mW/cm^3 reference level
PFfactor = 1;

NoMat = m1-1;
FreqFlag    = zeros(size(1:1:NoMat));
maxFreqPairs = floor((size(XCoreFreq,2))/2);  % or compute from your sheets
ConstantA      = zeros(NoMat, maxFreqPairs);
ConstantB      = zeros(NoMat, maxFreqPairs);
B_atPv_500     = zeros(NoMat, maxFreqPairs);
F_atPv_500     = zeros(NoMat, maxFreqPairs);
PF_atPv_500    = zeros(NoMat, maxFreqPairs);
beta_range     = zeros(NoMat, maxFreqPairs);
XCorePloss_3rd = zeros(NoMat, maxFreqPairs);
alpha_range    = zeros(NoMat, maxFreqPairs);
K1_range       = zeros(NoMat, maxFreqPairs);

for i = 1:1:NoMat
    DataSheetFreq = XCoreFreq(i,~isnan(XCoreFreq(i,:)));
    NoFreq = length(DataSheetFreq)/2;
    
    for j = 1:1:NoFreq
        % log10(P) = A*log10(B) + B0 at this data-sheet freq

        ConstantA(i,j) = (log10(XCorePloss(i,2*j)) - log10(XCorePloss(i,2*j-1)))/( ...
            log10(XCoreBfield(i,2*j)) - log10(XCoreBfield(i,2*j-1)));
        ConstantB(i,j) = log10(XCorePloss(i,2*j)) - ConstantA(i ,j)*log10(XCoreBfield( ...
            i,2*j));
        B_atPv_500(i,j) = 10^((log10(Pbar) - ConstantB(i,j))/ConstantA(i,j)); % in T
        F_atPv_500(i,j) = DataSheetFreq(2*j-1); % in Hz
        PF_atPv_500(i,j) = B_atPv_500(i,j)*F_atPv_500(i,j)^PFfactor;

        % For ./, it works, it just uses this to account for if frequency is a
        % sweeped variable i.e. an array.
        if (abs(fs_range - F_atPv_500(i,j))/fs_range <= 0.4)
            FreqFlag(i) = 1;
        end

        if (j > 1)
            % Steinmetz exponents across two adjacent frequencies
            beta_range(i,j) = log10(XCorePloss(i,2*j)/XCorePloss(i,2*j-1))/log10( ...
                XCoreBfield(i,2*j)/XCoreBfield(i,2*j-1));
            % Extrapolate third point on previous line (same B2)
            XCorePloss_3rd(i,j) = 10.^(ConstantA(i,j-1)*log10(XCoreBfield(i,2*j))+ ...
                ConstantB(i,j-1));
            alpha_range(i,j) = log10(XCorePloss_3rd(i,j)/XCorePloss(i,2*j))/log10( ...
                DataSheetFreq(2*j-3)/DataSheetFreq(2*j-1)); %(f2/f1)^alpha = P2/P1;
            K1_range(i,j) = XCorePloss(i,2*j)/(XCoreBfield(i,2*j)^beta_range(i,j))/( ...
                DataSheetFreq(2*j-1)^alpha_range(i,j)); %nW/cm3

            % Populate j-1 if missing (mirror j)
            if j == 2
                beta_range(i,j-1)  = beta_range(i,j);
                alpha_range(i,j-1) = alpha_range(i,j);
                K1_range(i,j-1) = XCorePloss(i,2*j-2)/(XCoreBfield(i,2*j-2)^ ...
                    beta_range(i,j-1))/(DataSheetFreq(2*j-3)^alpha_range(i,j-1));            
            end
        end
    end
end


%Core size
%-------------------------------------------

[m1,~] = size(raw);
TransformerCoreIndex = cell2mat(raw(2:m1,1));
XcoreVe = cell2mat(raw(2:m1,3))/(1000^3); % in m
XcoreAe = cell2mat(raw(2:m1,4))/(1000^2);
XcoreLe = cell2mat(raw(2:m1,5))/1000;

XcoreCoreShapeIndex = cell2mat(raw(2:m1,6));
XcorePriW = cell2mat(raw(2:m1,8))/1000;
XcorePriH = cell2mat(raw(2:m1,9))/1000;
XcoreSecW = cell2mat(raw(2:m1,10))/1000;
XcoreSecH = cell2mat(raw(2:m1,11))/1000;

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

%XcoreWindowH = 2*cell2mat(raw(2:m1,13))/1000;
%Replaced since GPT says the sheet already has full window height but this
%is based on my naming it so I have to cross reference with the thesis

ShuffleIndex = 1:1:length(TransformerCoreIndex);

% DESIGN SWEEP
%% ------------------------------------------------------------------------

CoreMatIndexSweep = find(FreqFlag);


% N-Dimensional grid expansion, 
% then shrinking based on some initial filters
% -------------------------------------

% Here I need to make this more efficient. Instad of making a 10
% dimensional grid to subsequently sweep, I should use a different method
% that is lighter on RAM.

[Po, fs, Vppeak, Vspeak, Vinsulation_max, matno_record, ShuffleXcoreIndex, Np, Mlp, Mls] = ndgrid( ...
    Po_range, fs_range, Vppeak_range, Vspeak_range, Vinsulation_max_range, ...
    CoreMatIndexSweep, ShuffleIndex, MinPriWinding:IncreNp:MaxPriWinding, ...
    1:IncreMlp:MaxMlp, 1:IncreMls:MaxMls);

% Flatten
Po = reshape(Po,[],1);
fs = reshape(fs,[],1);
Vppeak = reshape(Vppeak,[],1);
Vspeak = reshape(Vspeak,[],1);
Vinsulation_max = reshape(Vinsulation_max,[],1);
matno_record = reshape(matno_record,[],1);
Np = reshape(Np,[],1);
Mlp = reshape(Mlp,[],1);
Mls = reshape(Mls,[],1);
ui = XCoreMU(matno_record);
BSAT = XCoreBSAT(matno_record);
ShuffleXcoreIndex = reshape(ShuffleXcoreIndex,[],1);

% Map XcoreSize to actual size
Ve = XcoreVe(ShuffleXcoreIndex);
Ac = XcoreAe(ShuffleXcoreIndex);
W= XcoreWindowW(ShuffleXcoreIndex);
H= XcoreWindowH(ShuffleXcoreIndex);
Le = XcoreLe(ShuffleXcoreIndex);
PriW = XcorePriW(ShuffleXcoreIndex);
PriH = XcorePriH(ShuffleXcoreIndex);
SecW = XcoreSecW(ShuffleXcoreIndex);
SecH = XcoreSecH(ShuffleXcoreIndex);
XcoreIndex = TransformerCoreIndex(ShuffleXcoreIndex);
XcoreCoreShapeIndex = XcoreCoreShapeIndex(ShuffleXcoreIndex);

% Quick BSAT feasibility using datasheet formula
Bm_dummy = Vppeak/pi./fs./(2*Np.*Ac);
Keep_Bmindex = find(Bm_dummy < BSAT*BSAT_discount);
KeepIndex = Keep_Bmindex;

Po                   = Po(KeepIndex);
fs                   = fs(KeepIndex);
Vppeak               = Vppeak(KeepIndex);
Vspeak               = Vspeak(KeepIndex);
Vinsulation_max      = Vinsulation_max(KeepIndex);
matno_record         = matno_record(KeepIndex);
ui                   = ui(KeepIndex);
BSAT                 = BSAT(KeepIndex);
Ve                   = Ve(KeepIndex);
Ac                   = Ac(KeepIndex);
Le                   = Le(KeepIndex);
W                    = W(KeepIndex);
H                    = H(KeepIndex);
PriW                 = PriW(KeepIndex);
PriH                 = PriH(KeepIndex);
SecW                 = SecW(KeepIndex);
SecH                 = SecH(KeepIndex);
XcoreIndex           = XcoreIndex(KeepIndex);
XcoreCoreShapeIndex  = XcoreCoreShapeIndex(KeepIndex);
Np                   = Np(KeepIndex);
Mlp                  = Mlp(KeepIndex);
Mls                  = Mls(KeepIndex);


if isempty(Po)
    fprintf('Bm_dummy range: %.3f .. %.3f T\n', min(Bm_dummy), max(Bm_dummy));
    fprintf('BSAT*0.75 range: %.3f .. %.3f T\n', min(BSAT*BSAT_discount), max(BSAT*BSAT_discount));
    error(['Empty. Every candidate design violates the B-SAT screening.' ...
        'Increase transformer geometry size (likely), or Increase Np range, ' ...
        'or loosen other ranges, or ensure units match. Try again']);
    return; %#ok<UNRCH>
end


% Select Steinmetz params near each fs and expand design rows accordingly
FsnoNonzero = F_atPv_500(matno_record,:) > 0;
FsnoIndex = abs(fs - F_atPv_500(matno_record ,:))./fs <= 0.2;
% 20% away is the 0.2
matfsIndex = FsnoNonzero.*FsnoIndex;
matfs = F_atPv_500(matno_record,:).*matfsIndex;
K1 = K1_range(matno_record,:).*matfsIndex*1000; %convert from mW/cm3 to W/m3
alpha = alpha_range(matno_record,:).*matfsIndex;
beta = beta_range(matno_record,:).*matfsIndex;
[rowIdcs, ~] = find(matfs > 0);

% So far, each row of the above represent one DESIGN POINT (that has one
% set of electrical requirements, one core size, one core material, one Np, Mlp and Mls);

% Each row of matfs, Kl, alpha and beta also correspond to each DESIGN POINT
% However, they have more than one non-zero columns because each material
% may have more than one loss data points in their datasheets around the required frequency

% We need to expand the design point to incorporate different loss data
% points for one material.

% Find the indices of unique values in rowIdcs
[UniqueRowIdcs, ~] = unique(rowIdcs,'rows');
ColDuplicate = sum(matfs(UniqueRowIdcs,:)~=0,2);

% Test
if isempty(rowIdcs)
    error('CoreLoss:NoUsableColumns', ...
          ['No datasheet loss columns within ±40%% of fs=%g Hz for any selected material.\n' ...
           'Check CoreLossData.xlsx sheets: Freq/Bfield/Ploss pairs and fs_range.'], fs_range);
end

% Repeat by the number of loss data of each design point

Po                  = repelem(Po(UniqueRowIdcs), ColDuplicate);
fs                  = repelem(fs(UniqueRowIdcs), ColDuplicate);
Vppeak              = repelem(Vppeak(UniqueRowIdcs), ColDuplicate);
Vspeak              = repelem(Vspeak(UniqueRowIdcs), ColDuplicate);
Vinsulation_max     = repelem(Vinsulation_max(UniqueRowIdcs), ColDuplicate);
matno_record        = repelem(matno_record(UniqueRowIdcs), ColDuplicate);
ui                  = repelem(ui(UniqueRowIdcs), ColDuplicate);
BSAT                = repelem(BSAT(UniqueRowIdcs), ColDuplicate);
Ve                  = repelem(Ve(UniqueRowIdcs), ColDuplicate);
Ac                  = repelem(Ac(UniqueRowIdcs), ColDuplicate);
Le                  = repelem(Le(UniqueRowIdcs), ColDuplicate);
W                   = repelem(W(UniqueRowIdcs), ColDuplicate);
H                   = repelem(H(UniqueRowIdcs), ColDuplicate);
PriW                = repelem(PriW(UniqueRowIdcs), ColDuplicate);
PriH                = repelem(PriH(UniqueRowIdcs), ColDuplicate);
SecW                = repelem(SecW(UniqueRowIdcs), ColDuplicate);
SecH                = repelem(SecH(UniqueRowIdcs), ColDuplicate);
XcoreIndex          = repelem(XcoreIndex(UniqueRowIdcs), ColDuplicate);
XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex(UniqueRowIdcs), ColDuplicate);
Np                  = repelem(Np(UniqueRowIdcs), ColDuplicate);
Mlp                 = repelem(Mlp(UniqueRowIdcs), ColDuplicate);
Mls                 = repelem(Mls(UniqueRowIdcs), ColDuplicate);


% Collapse the nonzero columns into vectors aligned with repeats
matfs = nonzeros(reshape(matfs(UniqueRowIdcs,:)',[],1));
K1    = nonzeros(reshape(K1(UniqueRowIdcs,:)',   [],1));
beta  = nonzeros(reshape(beta(UniqueRowIdcs,:)', [],1));
alpha = nonzeros(reshape(alpha(UniqueRowIdcs,:)',[],1));

% For the remaining indexes, compute the more detailed parameters 
% like losses, size, subsequent weight
% -------------------------------------

if (isempty(Po))
    y = 0;
else

    % Here, a suitable wire size for primary and secondary is chosen to
    % meet current density limit and keep reasonable AC losses. That wire
    % size is then used to match everything else.

    % Electricals & losses
    % -------------------------------------
    skindepth = 1./sqrt(pi*fs*u0/rou);
    k  = Vspeak./Vppeak;
    Ns = round(Np.*k)+1;
    % Primary Current
    Iprms = Po/etaXfmer./(Vppeak/sqrt(2));
    Ippeak = Iprms*sqrt(2);
    % Secondary Current
    Isrms  = Po./(Vspeak/sqrt(2));
    Ispeak = Isrms*sqrt(2);
    % Required copper area

    % Flux and core loss
    % -------------------------------------

    lamda = Vppeak./pi./fs;
    % Equivalent load resistor
    Rload = Vspeak.*Vspeak./2./Po;
    % Input power (W)
    Pin = Po./etaXfmer;
    % Maximum total loss allowed (W)
    Ploss_est = Pin - Po;
    % Window area (m)
    Wa = H.*W;
    % Core volume (m3)
    Vcore = Ve;

    % Core weight (g)
    % -------------------------------------
    
    densityRow = CoreDensity(matno_record);
    assert(numel(densityRow)==numel(Vcore),'Core density mapping mismatch');
    Wcore = Vcore.*densityRow;
    
    % Calculate Bmax (T)
    % -------------------------------------
    Bm = lamda./(2.*Np.*Ac);
    
    % Calculate core loss (W)
    % -------------------------------------
    Pcore = CoreLossMultiple.*Vcore.*K1.*fs.^alpha.*Bm.^beta;

    % Determine wire type, and num of strands if Litz
    % -------------------------------------
    Areq_p=Iprms./Jwmax;                                % [m^2]
    Areq_s=Isrms./Jwmax;                                % [m^2]
    
    dsolid_p=2.*sqrt(Areq_p./pi);                       % [m]
    dsolid_s=2.*sqrt(Areq_s./pi);                       % [m]
    
    useSolid_p=(dsolid_p<=SolidGate.*(2.*skindepth));
    useSolid_s=(dsolid_s<=SolidGate.*(2.*skindepth));
    
    dstrand_litz=max(MinLitzDia,2.*skindepth);          % [m]
    Astrand=pi.*(dstrand_litz./2).^2;                   % [m^2]
    
    % Primary
    % ----------------
    Pri_Nstrands=ones(size(Iprms));
    Pri_Nstrands(~useSolid_p)=ceil(Areq_p(~useSolid_p)./Astrand(~useSolid_p));
    Pri_WireSize=dsolid_p;
    idx=~useSolid_p;
    if any(idx)
        Pri_WireSize(idx)=2.*sqrt((Pri_Nstrands(idx).*Astrand(idx))./(pi.*LitzFactor));
    end
    Pri_ds=dsolid_p;                                     % effective strand dia for AC model
    Pri_ds(idx)=dstrand_litz(idx);
    Pri_FullWireSize = Pri_WireSize + 2.*t_turn_p;
    
    % Secondary
    % ----------------
    Sec_Nstrands=ones(size(Isrms));
    Sec_Nstrands(~useSolid_s)=ceil(Areq_s(~useSolid_s)./Astrand(~useSolid_s));
    Sec_WireSize=dsolid_s;
    idx=~useSolid_s;
    if any(idx)
        Sec_WireSize(idx)=2.*sqrt((Sec_Nstrands(idx).*Astrand(idx))./(pi.*LitzFactor));
    end
    Sec_ds=dsolid_s;                                     % effective strand dia for AC model
    Sec_ds(idx)=dstrand_litz(idx);
    Sec_FullWireSize = Sec_WireSize + 2.*t_turn_s;

    A_pri_cu  = (pi.*(Pri_WireSize.^2))./4;      % m^2 per turn
    A_sec_cu  = (pi.*(Sec_WireSize.^2))./4;
    A_pri_full= (pi.*(Pri_FullWireSize.^2))./4;  % includes enamel
    A_sec_full= (pi.*(Sec_FullWireSize.^2))./4;
    
    CopperPacking  = (A_pri_cu.*Np + A_sec_cu.*(Ns./2))./ (H.*W);
    OverallPacking = (A_pri_full.*Np + A_sec_full.*(Ns./2))./ (H.*W);

    % Winding structure
    % ----------------
    
    DS_clearance = (UsePotting).*Potting_DS+(~UsePotting).*dielectricstrength_insulation;
    CoreInsulationThickness = CoronaMargin .* Vinsulation_max ./ DS_clearance;
    Pri_PerLayer=floor(Np./Mlp);
    Sec_PerLayer=floor(Ns./Mls);
    
    % Interlayer param.s
    % -----------------
    
    % Effective radial build added between adjacent secondary layers (geometry)
    t_interlayer_radial = Tape_Interlayer_Wraps .* Tape_Interlayer_Thickness ...
                    + (UsePotting .* Interlayer_Epoxy_Thickness);

    % Breakdown voltage capability of the tape + epoxy, derated
    Vbreak_interlayer = (Tape_Interlayer_Wraps .* Tape_Interlayer_Thickness .* Tape_Interlayer_DS) ...
                      + (UsePotting .* Interlayer_Epoxy_Thickness .* Potting_DS);
    Vbreak_interlayer = Vbreak_interlayer ./ CoronaMargin;

    % Estimated volts per layer
    V_per_layer_est = Vspeak ./ max(Mls,1);
    
    % Pass/fail index for interlayer stress
    Interlayer_index = find(V_per_layer_est <= Vbreak_interlayer);

    if useGPU
        % Move inputs for this block to GPU
        shape_g    = toGPU(int32(XcoreCoreShapeIndex));
        np_g       = toGPU(Np);
        ns_g       = toGPU(Ns);
        priW_g     = toGPU(PriW);
        priH_g     = toGPU(PriH);
        secW_g     = toGPU(SecW);
        secH_g     = toGPU(SecH);
        tcore_g    = toGPU(CoreInsulationThickness);
        mlp_g      = toGPU(Mlp);
        mls_g      = toGPU(Mls);
        priFull_g  = toGPU(Pri_FullWireSize);
        secFull_g  = toGPU(Sec_FullWireSize);
        wp_scalar  = toScalar(int32(Winding_Pattern));
        tinter_g = toGPU(t_interlayer_radial);
    
        % Element-wise compute for each candidate
        [TLp_g, TLs_g] = arrayfun(@computeTL_one,...
            shape_g, np_g, ns_g, priW_g, priH_g, secW_g, secH_g, ...
            tcore_g, mlp_g, mls_g, priFull_g, secFull_g, tinter_g, ...
            wp_scalar);
    
        TLp = toHost(TLp_g);
        TLs = toHost(TLs_g);
    else

        %CPU code version
      % Winding pattern of secondary
        isEE = (XcoreCoreShapeIndex == 1);
        isER = (XcoreCoreShapeIndex == 2);
        isU  = (XcoreCoreShapeIndex == 3);
        isUR = (XcoreCoreShapeIndex == 4);
    
        switch Winding_Pattern
            case 1 % -------- center leg --------
                % Secondary turns per layer
                Sec_PerLayer = floor(Ns./Mls);
    
                % Secondary group counts (per your table)
                Ns_group1 = floor((H - 2.*CoreInsulationThickness)./Sec_FullWireSize);
                Ns_group2 = zeros(size(Ns_group1));
                Ns_group3 = zeros(size(Ns_group1));
                Ns_group4 = zeros(size(Ns_group1));
                SupposeNs = Ns_group1 .* Mls;
    
                %% Total length of windings (primary unchanged)
                % EE center
                TLp(isEE) = Np(isEE).*2.*( ...
                    PriW(isEE) + PriH(isEE) + 4.*CoreInsulationThickness(isEE) + 2.*Mlp(isEE).*Pri_FullWireSize(isEE));
                TLs(isEE) = Ns(isEE).*2.*( ...
                    PriW(isEE) + PriH(isEE) + 4.*Mlp(isEE).*Pri_FullWireSize(isEE) + 8.*CoreInsulationThickness(isEE) ...
                    + 2.*( Mls(isEE).*Sec_FullWireSize(isEE) + max(Mls(isEE)-1,0).*t_interlayer_radial ) );
    
                % ER center
                TLp(isER) = 2.*pi.*Np(isER).*( ...
                    PriW(isER)./2 + CoreInsulationThickness(isER) + 0.5.*Mlp(isER).*Pri_FullWireSize(isER));
                TLs(isER) = 2.*pi.*Ns(isER).*( ...
                    PriW(isER)./2 + Mlp(isER).*Pri_FullWireSize(isER) + 2.*CoreInsulationThickness(isER) ...
                    + 0.5.*( Mls(isER).*Sec_FullWireSize(isER) + max(Mls(isER)-1,0).*t_interlayer_radial ) );
    
                % U center
                TLp(isU) = 2.*Np(isU).*( ...
                    PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) + 2.*Mlp(isU).*Pri_FullWireSize(isU));
                TLs(isU) = 2.*Ns(isU).*( ...
                    PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) ...
                    + 2.*( Mls(isU).*Sec_FullWireSize(isU) + max(Mls(isU)-1,0).*t_interlayer_radial ) );
    
                % UR center
                TLp(isUR) = 2.*pi.*Np(isUR).*( ...
                    PriW(isUR)./2 + CoreInsulationThickness(isUR) + 0.5.*Mlp(isUR).*Pri_FullWireSize(isUR));
                TLs(isUR) = 2.*pi.*Ns(isUR).*( ...
                    PriW(isUR)./2 + CoreInsulationThickness(isUR) ...
                    + 0.5.*( Mls(isUR).*Sec_FullWireSize(isUR) + max(Mls(isUR)-1,0).*t_interlayer_radial ) );
    
            case 2 % -------- double leg --------
                % Only consider half because symmetric
                Sec_PerLayer = floor(Ns./2./Mls);
    
                % Secondary group counts
                Ns_group1 = zeros(size(Sec_PerLayer));
                Ns_group2 = floor((W - 3.*CoreInsulationThickness - Mlp.*Pri_FullWireSize)./Sec_FullWireSize);
                Ns_group3 = floor((H - 2.*CoreInsulationThickness - 2.*Mls.*Sec_FullWireSize)./Sec_FullWireSize);
                Ns_group4 = Ns_group2;
                SupposeNs = Ns_group1.*Mls + Ns_group2.*Mls + Ns_group3.*Mls + Ns_group4.*Mls;
    
                %% Total length of windings
                % EE double
                TLp(isEE) = Np(isEE).*2.*( ...
                    PriW(isEE) + PriH(isEE) + 4.*CoreInsulationThickness(isEE) + 2.*Mlp(isEE).*Pri_FullWireSize(isEE));
                TLs(isEE) = Ns(isEE).*2.*( ...
                    SecW(isEE) + SecH(isEE) + 4.*CoreInsulationThickness(isEE) ...
                    + 2.*( Mls(isEE).*Sec_FullWireSize(isEE) + max(Mls(isEE)-1,0).*t_interlayer_radial ) );
    
                % ER double
                TLp(isER) = 2.*pi.*Np(isER).*( ...
                    PriW(isER)./2 + CoreInsulationThickness(isER) + 0.5.*Mlp(isER).*Pri_FullWireSize(isER));
                TLs(isER) = 2.*pi.*Ns(isER).*( ...
                    sqrt(2).*PriW(isER)./2 + CoreInsulationThickness(isER) ...
                    + 0.5.*( Mls(isER).*Sec_FullWireSize(isER) + max(Mls(isER)-1,0).*t_interlayer_radial ) );
    
                % U double
                TLp(isU) = 2.*Np(isU).*( ...
                    PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) + 2.*Mlp(isU).*Pri_FullWireSize(isU));
                TLs(isU) = 2.*Ns(isU).*( ...
                    PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) ...
                    + 2.*( Mls(isU).*Sec_FullWireSize(isU) + max(Mls(isU)-1,0).*t_interlayer_radial ) );
    
                % UR double
                TLp(isUR) = 2.*pi.*Np(isUR).*( ...
                    PriW(isUR)./2 + CoreInsulationThickness(isUR) + 0.5.*Mlp(isUR).*Pri_FullWireSize(isUR));
                TLs(isUR) = 2.*pi.*Ns(isUR).*( ...
                    PriW(isUR)./2 + CoreInsulationThickness(isUR) ...
                    + 0.5.*( Mls(isUR).*Sec_FullWireSize(isUR) + max(Mls(isUR)-1,0).*t_interlayer_radial ) );
    
            otherwise
                warning('Wrong winding pattern');
        end
    end



    % Calculate leakage inductance (not verified or used in this code)
    Lg = u0.*(W - Mls.*Sec_FullWireSize - Mlp.*Pri_FullWireSize).*SecH./H;
        %in Henry
    Xg = 2.*pi.*fs.*Lg;
    R_pri = Rload./Ns.^2;
    Lg_Lc_ratio = (W - 2.*CoreInsulationThickness - Mls.*Sec_FullWireSize - Mlp.* ...
        Pri_FullWireSize).*Le./ui./SecH./H;
    real_ratio = 1./(1 + Lg_Lc_ratio + Xg./R_pri);
    
    % Calculate Copper Loss
    %--------------------------------------------------------

    if useGPU
        % Move inputs to GPU
        TLp_g     = toGPU(TLp);
        TLs_g     = toGPU(TLs);
        PriDia_g  = toGPU(Pri_WireSize);
        SecDia_g  = toGPU(Sec_WireSize);
        PriN_g    = toGPU(Pri_Nstrands);
        SecN_g    = toGPU(Sec_Nstrands);
        dsPri_g   = toGPU(Pri_ds);   % effective strand dia (pri)
        dsSec_g   = toGPU(Sec_ds);   % effective strand dia (sec)
        skind_g   = toGPU(skindepth);
        Mlp_g     = toGPU(Mlp);
        Mls_g     = toGPU(Mls);
        rou_s     = toScalar(rou);   % scalar, stays on host
    
        % --- Dowell AC resistance ---
        Rdc_p_g = rou_s .* TLp_g ./ (pi .* (PriDia_g.^2) ./ 4);
        Rdc_s_g = rou_s .* TLs_g ./ (pi .* (SecDia_g.^2) ./ 4);
    
        % Primary
        Kp_g = sqrt(pi .* PriN_g) .* dsPri_g ./ (2 .* PriDia_g);
        xp_g = dsPri_g ./ (2 .* skind_g) .* sqrt(pi .* Kp_g);
    
        % Secondary
        Ks_g = sqrt(pi .* SecN_g) .* dsSec_g ./ (2 .* SecDia_g);
        xs_g = dsSec_g ./ (2 .* skind_g) .* sqrt(pi .* Ks_g);
    
        % Factors
        Fr_p_g = xp_g .* ((sinh(2 .* xp_g) + sin(2 .* xp_g)) ./ (cosh(2 .* xp_g) - cos(2 .* xp_g)) ...
                 + 2 .* ((Mlp_g.^2 .* PriN_g - 1) ./ 3) .* ((sinh(xp_g) - sin(xp_g)) ./ (cosh(xp_g) + cos(xp_g))));
        Fr_s_g = xs_g .* ((sinh(2 .* xs_g) + sin(2 .* xs_g)) ./ (cosh(2 .* xs_g) - cos(2 .* xs_g)) ...
                 + 2 .* ((Mls_g.^2 .* SecN_g - 1) ./ 3) .* ((sinh(xs_g) - sin(xs_g)) ./ (cosh(xs_g) + cos(xs_g))));
    
        Pri_Rac = gather(Rdc_p_g .* Fr_p_g);
        Sec_Rac = gather(Rdc_s_g .* Fr_s_g);
    
        % Copper loss
        Pcopper = (Iprms.^2 .* Pri_Rac + Isrms.^2 .* Sec_Rac);
    else
        PriKlayer = sqrt(pi.*Pri_Nstrands).*Pri_ds./2./(Pri_WireSize);
        Pri_xp = Pri_ds./2./skindepth.*sqrt(pi.*PriKlayer);
    
        SecKlayer = sqrt(pi.*Sec_Nstrands).*Sec_ds./2./(Sec_WireSize);
        Sec_xp = Sec_ds./2./skindepth.*sqrt(pi.*SecKlayer);
    
        TLp = TLp';
        TLs = TLs';
        Pri_Rdc = rou.*TLp./(pi.*Pri_WireSize.^2./4);
        Sec_Rdc = rou.*TLs./(pi.*Sec_WireSize.^2./4);
        Pri_Fr = Pri_xp.*((sinh(2.*Pri_xp) + sin(2.*Pri_xp))./(cosh(2.*Pri_xp) - cos(2.* ...
            Pri_xp)) + 2.*(Mlp.^2.*Pri_Nstrands - 1)./3.*(sinh(Pri_xp) - sin(Pri_xp))./( ...
            cosh(Pri_xp) + cos(Pri_xp)));
        Pri_Rac = Pri_Rdc.*Pri_Fr;
        Sec_Fr = Sec_xp.*((sinh(2.*Sec_xp) + sin(2.*Sec_xp))./(cosh(2.*Sec_xp) - cos(2.* ...
            Sec_xp)) + 2.*(Mls.^2.*Sec_Nstrands - 1)./3.*(sinh(Sec_xp) - sin(Sec_xp))./( ...
            cosh(Sec_xp) + cos(Sec_xp)));
        Sec_Rac = Sec_Rdc.*Sec_Fr;
        Pcopper = (Iprms.^2.*Pri_Rac + Isrms.^2.*Sec_Rac);
    end

    % Calculate temperature rise
    %----------------------------------------------------------------
    
    % Convert areas to cm^2 before applying the correlation
    Ac_cm2 = Ac .* 1e4;          % m^2 -> cm^2
    Wa_cm2 = Wa .* 1e4;          % m^2 -> cm^2
    Rth    = 16.31 .* (Ac_cm2 .* Wa_cm2) .^ (-0.405);   % K/W
    Tafterloss = Rth .* (Pcopper + Pcore) + 25;          % °C

    % Calculate the weight
    %----------------------------------------------------------------

    if useGPU
        shape_g   = toGPU(int32(XcoreCoreShapeIndex));
        H_g       = toGPU(H);
        W_g       = toGPU(W);
        priW_g    = toGPU(PriW);
        priH_g    = toGPU(PriH);
        tcore_g   = toGPU(CoreInsulationThickness);
        dens_g    = toGPU(CoreInsulationDensity);
    
        WeightCore_Insu = toHost(arrayfun(@coreInsu_one, ...
            shape_g, H_g, W_g, priW_g, priH_g, tcore_g, dens_g));
    else
        %CPU code version
        % Need to calculate WeightCore_Insu independently for each geometry
        % Preallocate
        WeightCore_Insu = zeros(size(Np));
        % EE: approximate insulating film area that wraps inner window walls
        WeightCore_Insu(isEE) = ( ...
            2.*H(isEE).*(PriW(isEE) + 2*PriH(isEE)) + ...   % vertical walls
            4.*W(isEE).*(PriW(isEE) + 2*PriH(isEE)) + ...   % horizontal walls
            H(isEE).*(2*PriW(isEE) + 2*PriH(isEE)) ...      % end caps / cheeks
            ).*CoreInsulationThickness(isEE).*CoreInsulationDensity;
    
        % ER: circular/rounded window (use circumference terms)
        WeightCore_Insu(isER) = ( ...
            sqrt(2)*pi.*H(isER).*PriW(isER) + ...           % vertical wrap
            sqrt(2)*pi.*2.*W(isER).*PriW(isER) + ...        % horizontal wrap
            H(isER).*pi.*PriW(isER) ...                     % cheeks
            ).*CoreInsulationThickness(isER).*CoreInsulationDensity;
    
        % U: rectangular, one side open; drop one horizontal wrap term
        WeightCore_Insu(isU) = ( ...
            2.*H(isU).*(PriW(isU) + 2*PriH(isU)) + ...
            2.*W(isU).*(PriW(isU) + 2*PriH(isU)) + ...
            H(isU).*(2*PriW(isU) + 2*PriH(isU)) ...
            ).*CoreInsulationThickness(isU).*CoreInsulationDensity;
    
        % UR: rounded U (use circular terms horizontally, rectangular vertically)
        WeightCore_Insu(isUR) = ( ...
            2.*H(isUR).*(PriW(isUR) + 2*PriH(isUR)) + ...
            pi.*W(isUR).*PriW(isUR) + ...                   % rounded horizontal
            H(isUR).*pi.*PriW(isUR) ...
            ).*CoreInsulationThickness(isUR).*CoreInsulationDensity;
    end
   
    % Wire weights
    WeightPri_copper = (pi.*Pri_WireSize.^2./4).*TLp.* CopperDensity;
    WeightPri_Insu = (pi.*(Pri_FullWireSize.^2 - Pri_WireSize.^2)./4).*TLp.*WireInsulationDensity;
    WeightSec_copper = (pi.*Sec_WireSize.^2./4).*TLs.*CopperDensity;
    WeightSec_Insu = (pi.*(Sec_FullWireSize.^2 - Sec_WireSize.^2)./4).*TLs.*WireInsulationDensity;

    % Compute total weight
    TotalWeight = Wcore + WeightPri_copper + WeightSec_copper + WeightPri_Insu + ...
        WeightSec_Insu + WeightCore_Insu;

    % Filter good designs
    %---------------------------------------------------------------
    
    % Filter by B
    B_index = find(Bm < BSAT*BSAT_discount);

    % Filter by Temperature
    P_loss_index = find(Pcopper + Pcore <= Ploss_est);
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);

    % Filter by weight
    TotalWeight_index = find(TotalWeight < MaxWeight);

    % Filter by packing factor
    OverallPackingmin_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactor);

    switch Winding_Pattern
        case 1 % center leg
            Ns_group1 = floor((H - 2*CoreInsulationThickness)./Sec_FullWireSize);
            Ns_group2 = zeros(size(Ns_group1));
            Ns_group3 = zeros(size(Ns_group1));
            Ns_group4 = zeros(size(Ns_group1));
            SupposeNs = Ns_group1.*Mls;
        case 2 % double leg
            Ns_group1 = zeros(size(Sec_PerLayer));
            Ns_group2 = floor((W - 3*CoreInsulationThickness - Mlp.*Pri_FullWireSize)./Sec_FullWireSize);
            Ns_group3 = floor((H - 2*CoreInsulationThickness - 2*Mls.*Sec_FullWireSize)./Sec_FullWireSize);
            Ns_group4 = Ns_group2;
            SupposeNs = Ns_group1.*Mls + Ns_group2.*Mls + Ns_group3.*Mls + Ns_group4.*Mls;
        otherwise
            error('Unsupported winding pattern');
    end

    if useGPU
        shape_g  = toGPU(int32(XcoreCoreShapeIndex));
        wp_s     = toScalar(int32(Winding_Pattern));
    
        PriPer_g = toGPU(Pri_PerLayer);
        PriFull_g= toGPU(Pri_FullWireSize);
        SecPer_g = toGPU(Sec_PerLayer);
        SecFull_g= toGPU(Sec_FullWireSize);
        Mlp_g    = toGPU(Mlp);
        Mls_g    = toGPU(Mls);
    
        Ns1_g    = toGPU(Ns_group1);
        Ns2_g    = toGPU(Ns_group2);
        Ns3_g    = toGPU(Ns_group3);
        Ns4_g    = toGPU(Ns_group4);
        H_g      = toGPU(H);
        W_g      = toGPU(W);
        Tcore_g  = toGPU(CoreInsulationThickness);
    
        SupposeNs_g= toGPU(SupposeNs);
        Ns_g      = toGPU(Ns);
    
        ok_g = arrayfun(@windowOk_one, shape_g, wp_s, ...
            PriPer_g, PriFull_g, SecPer_g, SecFull_g, ...
            Mlp_g, Mls_g, Ns1_g, Ns2_g, Ns3_g, Ns4_g, ...
            H_g, W_g, Tcore_g,tinter_g, SupposeNs_g, Ns_g);
    
        WindowFit_index = find(toHost(ok_g));

    else
        %CPU code version
        % Geometry masks
        is_EE_or_ER_core = (XcoreCoreShapeIndex==1) | (XcoreCoreShapeIndex==2);
        is_U_or_UR_core  = (XcoreCoreShapeIndex==3) | (XcoreCoreShapeIndex==4);
        
        % Winding pattern masks
        is_center_leg_pattern = (Winding_Pattern==1);
        is_double_leg_pattern = (Winding_Pattern==2);
        
        % Final pass/fail vector for this group of constraints
        meetsWindowFitConstraints = false(size(Np));
        
        % ===== EE/ER cores: center-leg pattern =====
        if is_center_leg_pattern
            % Height constraints (primary and secondary) and width constraint (both)
            EEER_center_pri_height_ok = is_EE_or_ER_core & ...
                (Pri_PerLayer.*Pri_FullWireSize <= H - 2.*CoreInsulationThickness);
        
            EEER_center_sec_height_ok = is_EE_or_ER_core & ...
                (Sec_PerLayer.*Sec_FullWireSize <= H - 2.*CoreInsulationThickness);
        
            EEER_center_width_ok = is_EE_or_ER_core & ...
                (Mlp.*Pri_FullWireSize + (Mls.*Sec_FullWireSize + max(Mls-1,0).*t_interlayer_radial) ...
                 <= W - 3.*CoreInsulationThickness);
        
            meetsWindowFitConstraints = meetsWindowFitConstraints | ...
                (EEER_center_pri_height_ok & EEER_center_sec_height_ok & EEER_center_width_ok);
        end
        
        % ===== EE/ER cores: double-leg pattern =====
        if is_double_leg_pattern
            EEER_double_pri_height_ok = is_EE_or_ER_core & ...
                (Pri_PerLayer.*Pri_FullWireSize <= H - 2.*CoreInsulationThickness);
        
            % Per table: (Ns3 + 2*Mls)*SecFullWireSize <= H - 2*TcoreInsu
            EEER_double_sec_height_ok = is_EE_or_ER_core & ...
                ((Ns_group3 + 2.*Mls).*Sec_FullWireSize <= H - 2.*CoreInsulationThickness);
        
            % Per table (width): Mlp*PriFullWireSize + Ns2*SecFullWireSize <= W - 3*TcoreInsu
            EEER_double_width_ok = is_EE_or_ER_core & ...
                (Mlp.*Pri_FullWireSize + Ns_group2.*Sec_FullWireSize <= W - 3.*CoreInsulationThickness);
        
            meetsWindowFitConstraints = meetsWindowFitConstraints | ...
                (EEER_double_pri_height_ok & EEER_double_sec_height_ok & EEER_double_width_ok);
        end
        
        % ===== U/UR cores (table gives one set; apply for either selected pattern) =====
        % Height: Pri_perlayer*PriFull + Ns1*SecFull <= H - 3*TcoreInsu
        UUR_height_ok = is_U_or_UR_core & ...
            (Pri_PerLayer.*Pri_FullWireSize + Ns_group1.*Sec_FullWireSize <= H - 3.*CoreInsulationThickness);
        
        % Width #1: Mlp*PriFull + Ns2*SecFull <= W - 3*TcoreInsu
        UUR_width1_ok = is_U_or_UR_core & ...
            (Mlp.*Pri_FullWireSize + Ns_group2.*Sec_FullWireSize <= W - 3.*CoreInsulationThickness);
        
        % Width #2: Mls*SecFull + Ns2*SecFull <= W - 2*TcoreInsu
        UUR_width2_ok = is_U_or_UR_core & ...
            ((Mls.*Sec_FullWireSize + max(Mls-1,0).*t_interlayer_radial) + Ns_group2.*Sec_FullWireSize ...
            <= W - 2.*CoreInsulationThickness);
        
        meetsWindowFitConstraints = meetsWindowFitConstraints | ...
            (UUR_height_ok & UUR_width1_ok & UUR_width2_ok);
        
        % Also require that the chosen layer/row plan can realize Ns overall
        hasEnoughSecondaryPositions = (SupposeNs >= Ns);
        meetsWindowFitConstraints = meetsWindowFitConstraints & hasEnoughSecondaryPositions;
        
        % Export indices used downstream
        WindowFit_index = find(meetsWindowFitConstraints);
    end

    % If no successful indexes, use step for every line here to see where
    % 0x1 array shows up, and thats the bottlenecking factor.
    Index_Meet_All = intersect(B_index,P_loss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tmin_index);
    Index_Meet_All = intersect(Index_Meet_All, Interlayer_index);
    Index_Meet_All = intersect(Index_Meet_All,TotalWeight_index);
    Index_Meet_All = intersect(Index_Meet_All,OverallPackingmin_index);
    Index_Meet_All = intersect(Index_Meet_All,OverallPackingmax_index);
    Index_Meet_All = intersect(Index_Meet_All, WindowFit_index);


    % Sort by total weight and keep only the lightest one
    [~,SortIndex] = sort(TotalWeight(Index_Meet_All));
    
    % Build Results Table 
    %% -----------------------------------------------------------------------------
    
    if(length(SortIndex) >= 1)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1:1));

        V_cu     = ((WeightPri_copper+WeightSec_copper))/CopperDensity;
        V_insu   = ((WeightCore_Insu+WeightPri_Insu+WeightSec_Insu))/CoreInsulationDensity;
        Volume_m3 = Vcore + V_cu + V_insu;

        Design(:, 1)  = Po(TotalWeightSortIndex);
        Design(:, 2)  = Vppeak(TotalWeightSortIndex);
        Design(:, 3)  = Vspeak(TotalWeightSortIndex);
        Design(:, 4)  = Vinsulation_max(TotalWeightSortIndex);
        Design(:, 5)  = fs(TotalWeightSortIndex);
        Design(:, 6)  = matno_record(TotalWeightSortIndex);
        Design(:, 7)  = matfs(TotalWeightSortIndex);
        Design(:, 8)  = Ac(TotalWeightSortIndex);
        Design(:, 9)  = H(TotalWeightSortIndex);
        Design(:,10)  = W(TotalWeightSortIndex);
        Design(:,11)  = Np(TotalWeightSortIndex);
        Design(:,12)  = Ns(TotalWeightSortIndex);
        Design(:,13)  = real_ratio(TotalWeightSortIndex);
        Design(:,14)  = Bm(TotalWeightSortIndex);
        Design(:,15)  = Pri_WireSize(TotalWeightSortIndex);
        Design(:,16)  = Pri_FullWireSize(TotalWeightSortIndex);
        Design(:,17)  = Sec_WireSize(TotalWeightSortIndex);
        Design(:,18)  = Sec_FullWireSize(TotalWeightSortIndex);
        Design(:,19) = Ippeak(TotalWeightSortIndex) ./ ...
                       (pi * Pri_Nstrands(TotalWeightSortIndex) .* Pri_ds(TotalWeightSortIndex).^2 / 4);
        Design(:,20) = Ispeak(TotalWeightSortIndex) ./ ...
                       (pi * Sec_Nstrands(TotalWeightSortIndex) .* Sec_ds(TotalWeightSortIndex).^2 / 4);
        Design(:,21) = Pri_Nstrands(TotalWeightSortIndex);
        Design(:,22) = Sec_Nstrands(TotalWeightSortIndex);
        Design(:,23) = Pri_PerLayer(TotalWeightSortIndex);
        Design(:,24) = Mlp(TotalWeightSortIndex);
        Design(:,25) = Sec_PerLayer(TotalWeightSortIndex);
        Design(:,26) = Mls(TotalWeightSortIndex);
        Design(:,27) = Ns_group1(TotalWeightSortIndex);
        Design(:,28) = Ns_group2(TotalWeightSortIndex);
        Design(:,29) = Ns_group3(TotalWeightSortIndex);
        Design(:,30) = Ns_group4(TotalWeightSortIndex);
        Design(:,31) = CopperPacking(TotalWeightSortIndex);
        Design(:,32) = OverallPacking(TotalWeightSortIndex);
        Design(:,33) = Pcore(TotalWeightSortIndex);
        Design(:,34) = Pcopper(TotalWeightSortIndex);
        Design(:,35) = Wcore(TotalWeightSortIndex);
        Design(:,36) = WeightPri_copper(TotalWeightSortIndex);
        Design(:,37) = WeightPri_Insu(TotalWeightSortIndex);
        Design(:,38) = WeightSec_copper(TotalWeightSortIndex);
        Design(:,39) = WeightSec_Insu(TotalWeightSortIndex);
        Design(:,40) = WeightCore_Insu(TotalWeightSortIndex);
        Design(:,41) = TotalWeight(TotalWeightSortIndex);
        Design(:,42) = Tafterloss(TotalWeightSortIndex);
        Design(:,43) = XcoreIndex(TotalWeightSortIndex);
        Design(:,44) = Volume_m3(TotalWeightSortIndex);
        y = Design;
    else
        y = zeros(1,44);
        disp('Requirements not met. Filtered indexes == 0');
    end
end
end



% ---------- TLp/TLs ----------
function [tLp,tLs] = computeTL_one(shape,np,ns,priW,priH,secW,secH,tcore,mlp,mls,priFull,secFull,tinter,wp)
% shape: 1=EE, 2=ER, 3=U, 4=UR; wp: 1=center-leg, 2=double-leg
% tinter: interlayer radial thickness per gap

    if wp==1
        % -------- Center-leg --------
        if shape==1   % EE rectangular window
            tLp = np*2*(priW+priH+4*tcore+2*mlp*priFull);
            % add 2*(mls*secFull + (mls-1)*tinter) for radial build
            tLs = ns*2*(priW+priH+4*mlp*priFull+8*tcore ...
                 + 2*( mls*secFull + max(mls-1,0)*tinter ));
        elseif shape==2 % ER circular
            tLp = 2*pi*np*(priW/2 + tcore + 0.5*mlp*priFull);
            tLs = 2*pi*ns*(priW/2 + mlp*priFull + 2*tcore ...
                 + 0.5*(mls*secFull + max(mls-1,0)*tinter));
        elseif shape==3 % U rectangular (center)
            tLp = 2*np*(priW+priH+4*tcore+2*mlp*priFull);
            tLs = 2*ns*(priW+priH+4*tcore ...
                 + 2*(mls*secFull + max(mls-1,0)*tinter));
        else            % UR circular (center)
            tLp = 2*pi*np*(priW/2 + tcore + 0.5*mlp*priFull);
            tLs = 2*pi*ns*(priW/2 + tcore + 0.5*(mls*secFull + max(mls-1,0)*tinter));
        end
    else
        % -------- Double-leg --------
        if shape==1   % EE rectangular (double)
            tLp = np*2*(priW+priH+4*tcore+2*mlp*priFull);
            tLs = ns*2*(secW+secH+4*tcore ...
                 + 2*(mls*secFull + max(mls-1,0)*tinter));
        elseif shape==2 % ER circular (double)
            tLp = 2*pi*np*(priW/2 + tcore + 0.5*mlp*priFull);
            tLs = 2*pi*ns*(sqrt(2)*priW/2 + tcore ...
                 + 0.5*(mls*secFull + max(mls-1,0)*tinter));
        elseif shape==3 % U rectangular (double)
            tLp = 2*np*(priW+priH+4*tcore+2*mlp*priFull);
            tLs = 2*ns*(priW+priH+4*tcore ...
                 + 2*(mls*secFull + max(mls-1,0)*tinter));
        else            % UR circular (double)
            tLp = 2*pi*np*(priW/2 + tcore + 0.5*mlp*priFull);
            tLs = 2*pi*ns*(priW/2 + tcore + 0.5*(mls*secFull + max(mls-1,0)*tinter));
        end
    end
end

% ---------- WeightCore_Insu ----------
function w = coreInsu_one(shape,H,W,priW,priH,tcore,dens)
    % area term depends on shape (your existing approximations)
    if shape==1          % EE: rectangular wraps both ways + cheeks
        area = 2*H*(priW+2*priH) + 4*W*(priW+2*priH) + H*(2*priW+2*priH);
    elseif shape==2      % ER: rounded
        area = sqrt(2)*pi*H*priW + sqrt(2)*pi*2*W*priW + H*pi*priW;
    elseif shape==3      % U: one side open; drop one horizontal wrap
        area = 2*H*(priW+2*priH) + 2*W*(priW+2*priH) + H*(2*priW+2*priH);
    else                 % UR: rounded horizontally, rectangular vertically
        area = 2*H*(priW+2*priH) + pi*W*priW + H*pi*priW;
    end
    w = area * tcore * dens;
end

% ---------- window-fit ----------
function ok = windowOk_one(shape,wp, PriPer, PriFull, SecPer, SecFull, ...
                           Mlp, Mls, Ns1, Ns2, Ns3, Ns4, H, W, Tcore, tinter,SupposeNs, Ns) %#ok<INUSD>
    if shape==1 || shape==2
        if wp==1
            % EE/ER center-leg
            okHpri = (PriPer*PriFull <= H - 2*Tcore);
            okHsec = (SecPer*SecFull <= H - 2*Tcore);
            % width includes Mls*SecFull plus (Mls-1)*tinter
            okW    = (Mlp*PriFull + (Mls*SecFull + max(Mls-1,0)*tinter) <= W - 3*Tcore);
            ok = okHpri && okHsec && okW;
        else
            % EE/ER double-leg (table width uses Ns2*SecFull only; no Mls term in that row)
            okHpri = (PriPer*PriFull <= H - 2*Tcore);
            okHsec = ((Ns3 + 2*Mls)*SecFull <= H - 2*Tcore);
            okW    = (Mlp*PriFull + Ns2*SecFull <= W - 3*Tcore);
            ok = okHpri && okHsec && okW;
        end
    else
        % U/UR (table gives two width constraints; only the one with Mls gets tinter)
        okH  = (PriPer*PriFull + Ns1*SecFull <= H - 3*Tcore);
        okW1 = (Mlp*PriFull + Ns2*SecFull   <= W - 3*Tcore);
        okW2 = ((Mls*SecFull + max(Mls-1,0)*tinter) + Ns2*SecFull <= W - 2*Tcore);
        ok = okH && okW1 && okW2;
    end
    ok = ok && (SupposeNs >= Ns);  % must be realizable
end

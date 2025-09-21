function y = Ecore_actual_EEER_xfmer_LCC(raw,raw1,raw2,raw3,raw4,raw5,raw6, ...
    Vppeak_range, Vspeak_range, Po_range, fs_range, Vinsulation_max_range, Winding_Pattern)


% Now I am going to edit this script so that it takes into account higher-voltage 
% transformer design parameters, such as volts-per-layer, epoxy potting, and interlayer tape.
% The inductor is much lower voltage since it's on the primary side, so it
% doesn't need these parameters.
UsePotting      = true;      % set true if vacuum-potted/encapsulated
Potting_DS      = 20e6;      % V/m, encapsulant dielectric strength (adjust per datasheet)
CoronaMargin    = 2.0;       % safety factor vs ideal breakdown
Vlayer_max      = 200;       % V allowed per secondary layer


%% Tunable Parameters

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
MaxPriWinding = 500;
% Incremental primary winding
IncreNp = 1;
% Maximum layer of primary winding
MaxMlp = 20;
% Incremental layer of primary winding
IncreMlp = 1;
% Maximum layer of secondary winding
MaxMls = 20;
% Incremental layer of secondary winding
IncreMls = 1;
% Minimum secondary wire diameter (m)
MinSecWireSize = 0.4e-3;             %#ok<NASGU>
% Max allowable transformer weight (g)
MaxWeight = 3000;
CoreInsulationDensity = 2.2e6;       % g/m^3 (Teflon)
WireInsulationDensity = 2.2e6;       % g/m^3 (Teflon)
% Saturation flux density derating
BSAT_discount = 0.9;
% Core loss multiplier
CoreLossMultiple = 1.5;
maxpackingfactor = 0.8;
minpackingfactor = 0.01;
% Litz copper fill
LitzFactor = 0.8;
% Not used directly here
BobbinWeightFactor = 0.5;             %#ok<NASGU>

%% Electrical constants

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

%% Body of function

% Save design results
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

%% ------------------------------------------------------------------------
% DESIGN SWEEP
CoreMatIndexSweep = find(FreqFlag);

% Here I need to make this more efficient. Instad of making a 10
% dimensional grid to subsequently sweep, I should use a different method
% that is lighter on RAM. Options include:
% - Preallocating RAM
% - 



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
KeepIndex = Keep_Bmindex; %#ok<NASGU>

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
ColDuplicate = sum(matfs(UniqueRowIdcs,:)~=0,2); %#ok<NASGU>

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

if (isempty(Po))
    y = 0;
else
    %Repeat elements by Primary Wire Number of Strands
    skindepth = 1./sqrt(pi*fs*u0/rou);
    %ds = max(skindepth ,MinLitzDia*ones(size(skindepth))); % take the skin depth litz
    ds = MinLitzDia*ones(size(skindepth));
    MinPriNstrands = floor((Po*2/etaXfmer./Vppeak/Jwmax)./(pi*ds.^2/4))+ 1;
    MaxPriNstrands = floor((Po*2/etaXfmer./Vppeak/Jwmax*1.0)./(pi*ds.^2/4))+ 1;
    

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
    Vppeak              = repelem(Vppeak, priCount);
    Vspeak              = repelem(Vspeak, priCount);
    Vinsulation_max     = repelem(Vinsulation_max, priCount);
    matno_record        = repelem(matno_record, priCount);
    ui                  = repelem(ui, priCount);
    BSAT                = repelem(BSAT, priCount);
    Ve                  = repelem(Ve, priCount);
    Ac                  = repelem(Ac, priCount);
    Le                  = repelem(Le, priCount);
    W                   = repelem(W, priCount);
    H                   = repelem(H, priCount);
    PriW                = repelem(PriW, priCount);
    PriH                = repelem(PriH, priCount);
    SecW                = repelem(SecW, priCount);
    SecH                = repelem(SecH, priCount);
    XcoreIndex          = repelem(XcoreIndex, priCount);
    XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex, priCount);
    Np                  = repelem(Np, priCount);
    Mlp                 = repelem(Mlp, priCount);
    Mls                 = repelem(Mls, priCount);
    matfs               = repelem(matfs, priCount);
    K1                  = repelem(K1, priCount);
    beta                = repelem(beta, priCount);
    alpha               = repelem(alpha, priCount);


    Pri_Nstrands = repmat((MinPriNstrands(1):1:MaxPriNstrands(1))',length( ...
        MaxPriNstrands),1);

    % Secondary: fixed to specific litz family (178-5790) -> 19 strands
    skindepth = 1./sqrt(pi*fs*u0/rou);
    ds = max(skindepth, MinLitzDia*ones(size(skindepth)));
    MinSecNstrands = 19*ones(size(skindepth));%178-5790 has 19 strands.
    MaxSecNstrands = 19*ones(size(skindepth));%178-5790 has 19 strands.
    secCount = MaxSecNstrands - MinSecNstrands + 1; %#ok<NASGU>
    
    Po                  = repelem(Po, secCount);
    fs                  = repelem(fs, secCount);
    Vppeak              = repelem(Vppeak, secCount);
    Vspeak              = repelem(Vspeak, secCount);
    Vinsulation_max     = repelem(Vinsulation_max, secCount);
    matno_record        = repelem(matno_record, secCount);
    ui                  = repelem(ui, secCount);
    BSAT                = repelem(BSAT, secCount);
    Ve                  = repelem(Ve, secCount);
    Ac                  = repelem(Ac, secCount);
    Le                  = repelem(Le, secCount);
    W                   = repelem(W, secCount);
    H                   = repelem(H, secCount);
    PriW                = repelem(PriW, secCount);
    PriH                = repelem(PriH, secCount);
    SecW                = repelem(SecW, secCount);
    SecH                = repelem(SecH, secCount);
    XcoreIndex          = repelem(XcoreIndex, secCount);
    XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex, secCount);
    Np                  = repelem(Np, secCount);
    Mlp                 = repelem(Mlp, secCount);
    Mls                 = repelem(Mls, secCount);
    matfs               = repelem(matfs, secCount);
    K1                  = repelem(K1, secCount);
    beta                = repelem(beta, secCount);
    alpha               = repelem(alpha, secCount);
    Pri_Nstrands        = repelem(Pri_Nstrands, secCount);

    Sec_Nstrands = repmat((MinSecNstrands(1):1:MaxSecNstrands(1))',length( ...
        MaxSecNstrands),1);


    % Electricals & losses
    k  = Vspeak./Vppeak;
    Ns = round(Np.*k)+1;
    % Primary Current
    Iprms = Po/etaXfmer./(Vppeak/sqrt(2));
    Ippeak = Iprms*sqrt(2);
    % Secondary Current
    Isrms  = Po./(Vspeak/sqrt(2));
    Ispeak = Isrms*sqrt(2);
    skindepth = 1./sqrt(pi*fs*u0/rou);
    ds = max(skindepth, MinLitzDia*ones(size(skindepth)));
    
    % Flux and core loss
    lamda = Vppeak./pi./fs;         % volt-second per radian
    % Check this value^

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
    Wcore = Vcore.*CoreDensity;
    % Calculate Bmax (T)
    Bm = lamda./(2.*Np.*Ac);
    % Calculate core loss (W)
    Pcore = CoreLossMultiple.*Vcore.*K1.*fs.^alpha.*Bm.^beta;
    
    % Wire
    % Primary wire diameter (m)
    Pri_WireSize = sqrt(Pri_Nstrands.*pi.*ds.^2./4./LitzFactor./pi).*2;
    % Primary wire diameter (m) including the insulation layer
    Pri_FullWireSize = Pri_WireSize + (Vppeak./dielectricstrength_insulation).*2;
    % Secondary wire diameter (m)
    Sec_WireSize = sqrt(Sec_Nstrands.*pi.*ds.^2./4./LitzFactor./pi).*2;
    % Secondary wire diameter (m) including the insulation layer
    Sec_FullWireSize = Sec_WireSize + (Vspeak./dielectricstrength_insulation/2).*2; %#ok<NASGU>
        % ^for above line, 100% dielectric strength, similar with Rubadue data
    % For 178-5790 only
    Sec_FullWireSize = 1.016./1000*ones(size(Sec_WireSize));
    CopperPacking = (pi.*Pri_WireSize.^2.*Np./4 + pi.*Sec_WireSize.^2.* ...
        Ns/2./4)./(H.*W);
    OverallPacking = (pi.*Pri_FullWireSize.^2.*Np./4 + pi.*Sec_FullWireSize ...
        .^2.*Ns/2./4)./(H.*W);

    % Winding structure
    CoreInsulationThickness  = Vinsulation_max./dielectricstrength_insulation;
    Pri_PerLayer=floor(Np./Mlp);
    Sec_PerLayer=floor(Ns./Mls);

    % Winding pattern of secondary
    % For ER, only allow center leg winding
    switch Winding_Pattern
        case 1 %center leg
            % Secondary turns per layer
            Sec_PerLayer = floor(Ns./Mls);
            % Number of rows in section 1, co-center with primary
            Ns_group1 = floor((H - 2*CoreInsulationThickness)./Sec_FullWireSize);
            Ns_group2 = zeros(size(Ns_group1));
            Ns_group3 = zeros(size (Ns_group1));
            Ns_group4 = zeros(size(Ns_group1));
            % Supposed secondary winding number if wind as mentioned above
            SupposeNs = Ns_group1.*Mls;
            %% Total length of windings
            TLp = Np.*2.*(PriW + PriH + 4*CoreInsulationThickness + 2.* ...
                Mlp.*Pri_FullWireSize);
            TLs= Ns.*2.*(PriW + PriH + 4*Mlp.*Pri_FullWireSize + 8* ...
                CoreInsulationThickness + 2.*Mls.*Sec_FullWireSize);
            % Recalculate XcoreCoreShapeIndex == 2, ER cores
            SelecIndex = find(XcoreCoreShapeIndex == 2);
            TLp(SelecIndex) = 2.*pi.*Np(SelecIndex).*(PriW(SelecIndex)./2 + ...
                CoreInsulationThickness(SelecIndex) + 0.5.*Mlp(SelecIndex) ...
                .*Pri_FullWireSize(SelecIndex));
            TLs(SelecIndex) = 2.*pi.*Ns(SelecIndex).*(PriW(SelecIndex)./2+ ...
            Mlp(SelecIndex).*Pri_FullWireSize(SelecIndex) + 2*CoreInsulationThickness ...
            (SelecIndex) + 0.5.*Mls(SelecIndex).*Sec_FullWireSize(SelecIndex));
        case 2 %double leg
            %% Winding pattern
            % Only consider half because symmetric
            Sec_PerLayer = floor(Ns./2./Mls);
            % Winding on double leg
            Ns_group1 = zeros(size(Sec_PerLayer)); % does not wind on center leg
            Ns_group2 = floor((W - 3*CoreInsulationThickness - Mlp.* ...
            Pri_FullWireSize)./Sec_FullWireSize);
            Ns_group3 = floor((H - 2*CoreInsulationThickness - 2*Mls.* ...
                Sec_FullWireSize)./Sec_FullWireSize);
            Ns_group4 = Ns_group2;
            SupposeNs = Ns_group1.*Mls + Ns_group2.*Mls + Ns_group3.*...
            Mls + Ns_group4.*Mls;
            %% Total length of windings
            TLp = Np.*2.*(PriW + PriH + 4*CoreInsulationThickness + 2.*Mlp.* ...
                Pri_FullWireSize);
            TLs = Ns.*2.*(SecW + SecH + 4*CoreInsulationThickness + 2.*Mls.* ...
                Sec_FullWireSize);
        otherwise
            warning('Wrong winding pattern');
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
    PriKlayer = sqrt(pi.*Pri_Nstrands).*ds./2./(Pri_WireSize);
    Pri_xp = ds./2./skindepth.*sqrt(pi.*PriKlayer);
    SecKlayer = sqrt(pi.*Sec_Nstrands).*ds./2./(Sec_WireSize);
    Sec_xp = ds./2./skindepth.*sqrt(pi.*SecKlayer);
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

    % Calculate temperature rise
    %----------------------------------------------------------------
    Rth = 16.31.*1e-3.*(Ac.*Wa).^(-0.405);
    
    
    %Rth = 16.31.*((Ac*1e4).*(Wa*1e4)).^(-0.405);
    %(GPT says this is a unit mismatch and to replace with above,
    % since 16.31 is in cm, might be unit mismatch)

    Tafterloss = Rth.*(Pcopper + Pcore) + 25;

    % Calculate the weight
    %----------------------------------------------------------------
    WeightPri_copper = pi.*Pri_WireSize.^2./4.*TLp.* CopperDensity;
    WeightPri_Insu = pi.*(Pri_FullWireSize.^2 - Pri_WireSize.^2)./4.*TLp.* ...
        WireInsulationDensity;
    WeightSec_copper = pi.*Sec_WireSize.^2./4.*TLs.*CopperDensity;
    WeightSec_Insu = pi.*(Sec_FullWireSize.^2 - Sec_WireSize.^2)./4.*TLs.* ...
        WireInsulationDensity;
    WeightCore_Insu =(2.*H.*(PriW + 2*PriH) + 4.*W.*(PriW + 2*PriH) + H.*(2*PriW + ...
        2*PriH)).*CoreInsulationThickness.*CoreInsulationDensity;

    % Recalculate XcoreCoreShapeIndex == 2, ER cores
    SelecIndex = find(XcoreCoreShapeIndex == 2);
    WeightCore_Insu(SelecIndex) = (sqrt(2)*pi*H(SelecIndex).*PriW(SelecIndex) + ...
        sqrt(2)*pi*2*W(SelecIndex).*PriW(SelecIndex) + H(SelecIndex)*pi.*PriW( ...
            SelecIndex)).*CoreInsulationThickness(SelecIndex).*CoreInsulationDensity;
    TotalWeight = Wcore + WeightPri_copper + WeightSec_copper + WeightPri_Insu + ...
        WeightSec_Insu + WeightCore_Insu;

    % Filter good designs
    %---------------------------------------------------------------
    B_index = find(Bm < BSAT*BSAT_discount);
    P_loss_index = find(Pcopper + Pcore <= Ploss_est);
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);
    TotalWeight_index = find(TotalWeight < MaxWeight);

    OverallPackingmin_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactor);

    Mlp_index = find(Mlp.*Pri_FullWireSize <=W - 3*CoreInsulationThickness);
    Pri_PerLayer_index = find(Pri_PerLayer.*Pri_FullWireSize < H - 2* ...
        CoreInsulationThickness);

    % make sure pri and sec has enough insulation in between
    switch Winding_Pattern
        case 1 %center leg
            Mls_index = find(Mls.*Sec_FullWireSize + Mlp.*Pri_FullWireSize ...
                <=W - 3*CoreInsulationThickness);
        case 2 % double leg, Ns2 or Ns4 together with primary fits in the window 
            % width, Ns3 and primary also fits in the window width
            Mls_index = intersect(find(Ns_group2.*Sec_FullWireSize + Mlp.* ...
                Pri_FullWireSize <= W - 3*CoreInsulationThickness),...
                find(Ns_group3.*Sec_FullWireSize + Mlp.*Pri_FullWireSize ...
                <=W - 3*CoreInsulationThickness));
        otherwise
            disp('Winding pattern does not meet Mls requirement');
    end


    % Secondary winding index
    Ns_group1_index = find(Ns_group1 >= 0);
    Ns_group2_index = find(Ns_group2 >= 0);
    Ns_group3_index = find(Ns_group3 >= 0);
    Ns_group4_index = find(Ns_group4 >= 0);
    SupposeNs_index = find(SupposeNs >= Ns);

    Index_Meet_All = intersect(B_index,P_loss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tmin_index);
    Index_Meet_All = intersect(Index_Meet_All,TotalWeight_index);
    Index_Meet_All = intersect(Index_Meet_All,OverallPackingmin_index);
    Index_Meet_All = intersect(Index_Meet_All,OverallPackingmax_index);
    Index_Meet_All = intersect(Index_Meet_All,Mlp_index);
    Index_Meet_All = intersect(Index_Meet_All,Pri_PerLayer_index);
    Index_Meet_All = intersect(Index_Meet_All,Mls_index);
    Index_Meet_All = intersect(Index_Meet_All,Ns_group1_index);
    Index_Meet_All = intersect(Index_Meet_All,Ns_group2_index);
    Index_Meet_All = intersect(Index_Meet_All,Ns_group3_index);
    Index_Meet_All = intersect(Index_Meet_All,Ns_group4_index);
    Index_Meet_All = intersect(Index_Meet_All,SupposeNs_index);

    % Sort by total weight and keep only the lightest one
    %
    [~,SortIndex] = sort(TotalWeight(Index_Meet_All));

    if(length(SortIndex) >= 1)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1:1));

        V_core   = (Wcore)/CoreDensity;
        V_cu     = ((WeightPri_copper+WeightSec_copper))/CopperDensity;
        V_insu   = ((WeightCore_Insu+WeightPri_Insu+WeightSec_Insu))/CoreInsulationDensity;
        Volume_m3 = V_core + V_cu + V_insu;

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
                       (pi * Pri_Nstrands(TotalWeightSortIndex) .* ds(TotalWeightSortIndex).^2 / 4);
        Design(:,20) = Ispeak(TotalWeightSortIndex) ./ ...
                       (pi * Sec_Nstrands(TotalWeightSortIndex) .* ds(TotalWeightSortIndex).^2 / 4);
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

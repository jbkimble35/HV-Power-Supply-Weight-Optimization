function y = Ecore_actual_EEER_xfmer_LCC(raw,raw1,raw2,raw3,raw4,raw5,raw6, ...
    Vppeak_range, Vspeak_range, Po_range, fs_range, Vinsulation_max_range,...
    Winding_Pattern,layerTapeUse,enamelThickness,kaptonDielStrength,kaptonThickness,...
    MinTapeMargin,kaptonDensity,CoreInsulationDensity,WireInsulationDensity,dielectricstrength_insulation, ...
    etaXfmer,TmaxX,TminX,MinPriWindingX,MaxPriWindingX,IncreNpX,MaxMlpX,IncreMlpX,MaxMlsX,IncreMlsX,MaxWeightX, ...
    BSAT_discountX,CoreLossMultipleX,maxpackingfactorX,minpackingfactorX,LitzFactor,MinWireDia,...
    Jwmax,MinLitzDia,CopperDensity,rou,u0)

% Body of function
%% --------------------------------------------------------------------------------------

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
XCoreInitialMU    = cell2mat(raw5(2:m1,3));
[m1,~]    = size(raw6);
CoreDensity   = cell2mat(raw6(2:m1,3))*1000000;

% Build Steinmetz parameter sets around the target frequency
% ----------------------------------------------

Pbar = 500;      % mW/cm^3 reference level
PFfactor = 1;
NoMat = m1-1;
FreqFlag    = zeros(size(1:1:NoMat));
maxFreqPairs = floor((size(XCoreFreq,2))/2);
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

        if (abs(fs_range - F_atPv_500(i,j))/fs_range <= 0.2)
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


% Core size
%-------------------------------------------

% CoreSizeData should all be in mm, mm2, and mm3.
[m1,~] = size(raw);
TransformerCoreIndex = cell2mat(raw(2:m1,1));
XcoreVe = cell2mat(raw(2:m1,3))/(1000^3); % mm3 to m3
XcoreAe = cell2mat(raw(2:m1,4))/(1000^2); % mm2 to m2
XcoreLe = cell2mat(raw(2:m1,5))/1000; % mm to m

XcoreCoreShapeIndex = cell2mat(raw(2:m1,6));
XcorePriW = cell2mat(raw(2:m1,8))/1000; % mm to m
XcorePriH = cell2mat(raw(2:m1,9))/1000; % mm to m
XcoreSecW = cell2mat(raw(2:m1,10))/1000; % mm to m
XcoreSecH = cell2mat(raw(2:m1,11))/1000; % mm to m

XcoreWindowW = cell2mat(raw(2:m1,12))/1000;
XcoreWindowH = cell2mat(raw(2:m1,13))/1000;
ShuffleIndex = 1:1:length(TransformerCoreIndex);

% DESIGN SWEEP
%% ------------------------------------------------------------------------

CoreMatIndexSweep = find(FreqFlag);

% N-Dimensional grid expansion, 
% then shrinking based on some initial filters
% -------------------------------------

[Po, fs, Vppeak, Vspeak, Vinsulation_max, matno_record, ShuffleXcoreIndex, Np, Mlp, Mls] = ndgrid( ...
    Po_range, fs_range, Vppeak_range, Vspeak_range, Vinsulation_max_range, ...
    CoreMatIndexSweep, ShuffleIndex, MinPriWindingX:IncreNpX:MaxPriWindingX, ...
    1:IncreMlpX:MaxMlpX, 1:IncreMlsX:MaxMlsX);

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
ui = XCoreInitialMU(matno_record);
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
Bm_dummy = Vppeak./(pi.*fs.*2.*Np.*Ac);
Keep_Bmindex = find(Bm_dummy < BSAT*BSAT_discountX);
KeepIndex = Keep_Bmindex;

% 1st Filter
% Applying initial BSAT filter
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
    fprintf('BSAT*0.75 range: %.3f .. %.3f T\n', min(BSAT*BSAT_discountX), max(BSAT*BSAT_discountX));
    error(['Empty. Every candidate design violates the B-SAT screening.' ...
        'Increase transformer geometry size (likely), or Increase Np range, ' ...
        'or loosen other ranges, or ensure units match. Try again']);
    return; %#ok<UNRCH>
end


% Select Steinmetz params near each fs and expand design rows accordingly
FsnoNonzero = F_atPv_500(matno_record,:) > 0;
FsnoIndex = abs(fs - F_atPv_500(matno_record ,:))./fs <= 0.4; % 40% away is the 0.4
matfsIndex = FsnoNonzero.*FsnoIndex;
matfs = F_atPv_500(matno_record,:).*matfsIndex;
K1 = K1_range(matno_record,:).*matfsIndex*1000; %convert from mW/cm3 to W/m3
alpha = alpha_range(matno_record,:).*matfsIndex;
beta = beta_range(matno_record,:).*matfsIndex;
[rowIdcs, ~] = find(matfs > 0);

% So far, each row of the above represent one DESIGN POINT (that has one
% set of electrical requirements, one core size, one core material, one Np, Mlp and Mls);

% Each row of matfs, K1, alpha and beta also correspond to each DESIGN POINT
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
          ['No datasheet loss columns within ±20% of fs=%g Hz for any selected material.\n' ...
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
% -----------------------------------------------------------------

if (isempty(Po))
    y = 0;
else

    % Here, a suitable wire size for primary and secondary is chosen to
    % meet current density limit and keep reasonable AC losses. That wire
    % size is then used to match everything else.

    % The wire size varies from the original thesis script as litz wire is
    % not automatically used. Since the thesis script intended for use with
    % higher frequency cores, litz was needed, but here the number of
    % strands necessary for skindepth to be throughout the wire is just computed

    % Electricals & losses
    % -------------------------------------

    skindepth = 1./sqrt(pi*fs*u0/rou);
    Ns = ceil(Np.*(Vspeak./Vppeak));
    % Primary Current
    Iprms = (Po/etaXfmer)./(Vppeak/sqrt(2));
    Ippeak = Iprms*sqrt(2);
    % Secondary Current
    Isrms  = Po./(Vspeak/sqrt(2));
    Ispeak = Isrms*sqrt(2);
    % Required copper area

    % Core weight (g)
    % -------------------------------------
    
    Wcore = Ve.*CoreDensity(matno_record);
    
    % Calculate Bmax (T)
    % -------------------------------------

    Bm = (Vppeak./pi./fs)./(2.*Np.*Ac);
    
    % Calculate core loss (W)
    % -------------------------------------

    Pcore = CoreLossMultipleX.*Ve.*K1.*fs.^alpha.*Bm.^beta;

    % Determine wire type, and num of strands if Litz
    % -----------------------------------------------

    Areq_p=Iprms./Jwmax;                                % [m^2]
    Areq_s=Isrms./Jwmax;                                % [m^2]
    
    dsolid_p=2.*sqrt(Areq_p./pi);                       % [m]
    dsolid_s=max(MinWireDia,2.*sqrt(Areq_s./pi));                       % [m]
    
    useSolid_p=(dsolid_p<=2.*skindepth);
    useSolid_s=(dsolid_s<=2.*skindepth);
    
    % Computes for skindepth but not DC resistance...
    MinLitzDia = MinLitzDia*ones(length(skindepth),1);
    dstrandMinlitz=max(MinLitzDia,2.*skindepth);          % [m]
    AstrandMin=pi.*(dstrandMinlitz./2).^2;                   % [m^2]
    
    % Primary wire diameter, number of strands, and jacket+wire diameter
    % ----------------

    Pri_Nstrands=ones(size(Iprms));
    Pri_Nstrands(~useSolid_p)=ceil(Areq_p(~useSolid_p)./AstrandMin(~useSolid_p));

    Pri_WireDia=max(MinWireDia,dsolid_p);
    idx=~useSolid_p;
    if any(idx)
        Pri_WireDia(idx)=2.*sqrt(Pri_Nstrands(idx).*AstrandMin(idx).*LitzFactor./pi);
    end
    % Strand diameter
    Pri_ds=max(MinWireDia,dsolid_p);
    Pri_ds(idx)=dstrandMinlitz(idx);

    Pri_FullWireDia = Pri_WireDia + 2.*Vppeak./dielectricstrength_insulation;
    if layerTapeUse
        Pri_FullWireDia = Pri_WireDia + enamelThickness.*2;
    end

    % Secondary wire diameter, number of strands, and jacket+wire diameter
    % ----------------

    Sec_Nstrands=ones(size(Isrms));
    Sec_Nstrands(~useSolid_s)=ceil(Areq_s(~useSolid_s)./AstrandMin(~useSolid_s));
    Sec_WireDia=max(MinWireDia,dsolid_s);
    idx=~useSolid_s;
    if any(idx)
        Sec_WireDia(idx)=2.*sqrt(Sec_Nstrands(idx).*AstrandMin(idx).*LitzFactor./pi);
    end
    % Strand diameter
    Sec_ds=max(MinWireDia,dsolid_s);
    Sec_ds(idx)=dstrandMinlitz(idx);

    % Changed vsp/die to have 2* factor like pri, then reverted.
    Sec_FullWireDia = Sec_WireDia + Vspeak./dielectricstrength_insulation;
    if layerTapeUse
        Sec_FullWireDia = Sec_WireDia + enamelThickness.*2;
    end

    % Winding structure
    % --------------------
    
    CoreInsulationThickness = Vinsulation_max./dielectricstrength_insulation;
    % Why is this calculated here and in the winding pattern part? I
    % deleted that one, which also had Sec_PerLayer calculated as half
    Pri_PerLayer=floor(Np./Mlp);
    Sec_PerLayer=floor(Ns./Mls);

    if layerTapeUse
        numTapePerLayerPri = ceil((Vppeak./Mlp)./(kaptonDielStrength*kaptonThickness));
        numTapePerLayerSec = ceil((Vinsulation_max./Mls)./(kaptonDielStrength*kaptonThickness));
        tTapePri = max(Mlp-1,0).*numTapePerLayerPri.*kaptonThickness;
        tTapeSec = max(Mls-1,0).*numTapePerLayerSec.*kaptonThickness;
    else
        numTapePerLayerPri = zeros(size(Mlp));
        numTapePerLayerSec = zeros(size(Mls));
        tTapePri = zeros(size(Mlp));
        tTapeSec = zeros(size(Mls));
    end
    % Computes mean length of turn for pri and sec, accounting for geometry
    % and winding pattern
    %-----------------------------------------------------------------------------
    
    % Winding pattern of secondary
    isEE = (XcoreCoreShapeIndex == 1);
    isER = (XcoreCoreShapeIndex == 2);
    isU  = (XcoreCoreShapeIndex == 3);
    isUR = (XcoreCoreShapeIndex == 4);

    switch Winding_Pattern
        case 1 % center leg ----------------

            Ns_group1 = floor((H - 2.*CoreInsulationThickness)./Sec_FullWireDia);
            Ns_group2 = zeros(size(Ns_group1));
            Ns_group3 = zeros(size(Ns_group1));
            Ns_group4 = zeros(size(Ns_group1));
            SupposeNs = Ns_group1 .* Mls;

            %% Total length of windings (primary unchanged)
            % EE center
            TLp(isEE) = Np(isEE).*2.*( ...
                PriW(isEE) + PriH(isEE) + 4.*CoreInsulationThickness(isEE) ...
                + 2.*Mlp(isEE).*Pri_FullWireDia(isEE) + 2.*tTapePri(isEE) );
            TLs(isEE) = Ns(isEE).*2.*( ...
                PriW(isEE) + PriH(isEE) ...
                + 4.*Mlp(isEE).*Pri_FullWireDia(isEE) + 4.*tTapePri(isEE) ...
                + 8.*CoreInsulationThickness(isEE) ...
                + 2.*( Mls(isEE).*Sec_FullWireDia(isEE) + tTapeSec(isEE) ) );

            % ER center
            TLp(isER) = 2.*pi.*Np(isER).*( ...
                PriW(isER)./2 + CoreInsulationThickness(isER) ...
                + 0.5.*Mlp(isER).*Pri_FullWireDia(isER) + 0.5.*tTapePri(isER) );
            TLs(isER) = 2.*pi.*Ns(isER).*( ...
                PriW(isER)./2 ...
                + Mlp(isER).*Pri_FullWireDia(isER) + tTapePri(isER) ...
                + 2.*CoreInsulationThickness(isER) ...
                + 0.5.*( Mls(isER).*Sec_FullWireDia(isER) + tTapeSec(isER) ) );

            % U center
            TLp(isU) = 2.*Np(isU).*( ...
                PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) ...
                + 2.*Mlp(isU).*Pri_FullWireDia(isU) + 2.*tTapePri(isU) );
            TLs(isU) = 2.*Ns(isU).*( ...
                PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) ...
                + 2.*( Mls(isU).*Sec_FullWireDia(isU) + tTapeSec(isU) ) );

            % UR center
            TLp(isUR) = 2.*pi.*Np(isUR).*( ...
                PriW(isUR)./2 + CoreInsulationThickness(isUR) ...
                + 0.5.*Mlp(isUR).*Pri_FullWireDia(isUR) + 0.5.*tTapePri(isUR) );
            TLs(isUR) = 2.*pi.*Ns(isUR).*( ...
                PriW(isUR)./2 + CoreInsulationThickness(isUR) ...
                + 0.5.*( Mls(isUR).*Sec_FullWireDia(isUR) + tTapeSec(isUR) ) );

        case 2 % double leg ----------------

            % Secondary group counts
            Ns_group1 = zeros(size(Sec_PerLayer));
            Ns_group2 = floor((W - 3.*CoreInsulationThickness - Mlp.*Pri_FullWireDia)./Sec_FullWireDia);
            Ns_group3 = floor((H - 2.*CoreInsulationThickness - 2.*Mls.*Sec_FullWireDia)./Sec_FullWireDia);
            Ns_group4 = Ns_group2;
            SupposeNs = Ns_group1.*Mls + Ns_group2.*Mls + Ns_group3.*Mls + Ns_group4.*Mls;

            %% Total length of windings
            % EE double
            TLp(isEE) = Np(isEE).*2.*( ...
                PriW(isEE) + PriH(isEE) + 4.*CoreInsulationThickness(isEE) ...
                + 2.*Mlp(isEE).*Pri_FullWireDia(isEE) + 2.*tTapePri(isEE) );
            TLs(isEE) = Ns(isEE).*2.*( ...
                SecW(isEE) + SecH(isEE) + 4.*CoreInsulationThickness(isEE) ...
                + 2.*( Mls(isEE).*Sec_FullWireDia(isEE) + tTapeSec(isEE) ) );

            % ER double
            TLp(isER) = 2.*pi.*Np(isER).*( ...
                PriW(isER)./2 + CoreInsulationThickness(isER) ...
                + 0.5.*Mlp(isER).*Pri_FullWireDia(isER) + 0.5.*tTapePri(isER) );
            TLs(isER) = 2.*pi.*Ns(isER).*( ...
                sqrt(2).*PriW(isER)./2 + CoreInsulationThickness(isER) ...
                + 0.5.*( Mls(isER).*Sec_FullWireDia(isER) + tTapeSec(isER) ) );

            % U double
            TLp(isU) = 2.*Np(isU).*( ...
                PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) ...
                + 2.*Mlp(isU).*Pri_FullWireDia(isU) + 2.*tTapePri(isU) );
            TLs(isU) = 2.*Ns(isU).*( ...
                PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU) ...
                + 2.*( Mls(isU).*Sec_FullWireDia(isU) + tTapeSec(isU) ) );

            % UR double
            TLp(isUR) = 2.*pi.*Np(isUR).*( ...
                PriW(isUR)./2 + CoreInsulationThickness(isUR) ...
                + 0.5.*Mlp(isUR).*Pri_FullWireDia(isUR) + 0.5.*tTapePri(isUR) );
            TLs(isUR) = 2.*pi.*Ns(isUR).*( ...
                PriW(isUR)./2 + CoreInsulationThickness(isUR) ...
                + 0.5.*( Mls(isUR).*Sec_FullWireDia(isUR) + tTapeSec(isUR) ) );

        otherwise
            warning('Wrong winding pattern');
    end

    % Calculate Copper Loss
    %--------------------------------------------------------

    PriKlayer = sqrt(pi.*Pri_Nstrands).*Pri_ds./(2.*Pri_WireDia);
    Pri_xp = Pri_ds./(sqrt(pi.*PriKlayer).*2.*skindepth);

    SecKlayer = sqrt(pi.*Sec_Nstrands).*Sec_ds./(2.*Sec_WireDia);
    Sec_xp = Sec_ds./(sqrt(pi.*SecKlayer).*2.*skindepth);

    TLp = TLp';
    TLs = TLs';

    % Fixed overestimation by litzfactor
    Pri_Rdc = rou .* TLp ./ (pi.*dsolid_p.^2./4);
    Sec_Rdc = rou .* TLs ./ (pi.*dsolid_s.^2./4);
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
    Wa=2.*H.*W.*(isEE|isER)+H.*W.*(isU|isUR);
    Rth    = 16.31e-3.*(Ac.*Wa).^(-0.405);   % K/W
    Tafterloss = Rth.*(Pcopper+Pcore)+25;          % °C

    % Weight of core insulation
    %----------------------------------------------------------------

    % This is just core liner of teflon

    WeightCore_Insu = zeros(size(Np));
    % EE
    WeightCore_Insu(isEE) = ( ...
        2.*H(isEE).*(PriW(isEE) + 2*PriH(isEE)) + ...  
        4.*W(isEE).*(PriW(isEE) + 2*PriH(isEE)) + ...   
        H(isEE).*(2*PriW(isEE) + 2*PriH(isEE)) ...      
        ).*CoreInsulationThickness(isEE).*CoreInsulationDensity;

    % ER
    WeightCore_Insu(isER) = ( ...
        sqrt(2)*pi.*H(isER).*PriW(isER) + ...          
        sqrt(2)*pi.*2.*W(isER).*PriW(isER) + ...        
        H(isER).*pi.*PriW(isER) ...                     
        ).*CoreInsulationThickness(isER).*CoreInsulationDensity;

    % U
    WeightCore_Insu(isU) = ( ...
        2.*H(isU).*(2.*PriW(isU) + 2.*PriH(isU)) + ...
        2.*W(isU).*(2.*PriW(isU) + 2.*PriH(isU)) ) ...
        .* CoreInsulationThickness(isU) .* CoreInsulationDensity;


    % UR
    WeightCore_Insu(isUR) = ( ...
        2*H(isUR).*(pi.*PriW(isUR)) + 2*W(isUR).*(pi.*PriW(isUR))).*...
        CoreInsulationThickness(isUR).*CoreInsulationDensity;

   
    % Wire weights
    %---------------------------------------------

    WeightPri_copper = pi.*Pri_WireDia.^2./4.*TLp.*CopperDensity;
    WeightPri_Insu = (pi.*(Pri_FullWireDia.^2 - Pri_WireDia.^2)./4).*TLp.*WireInsulationDensity;
    WeightSec_copper = pi.*Sec_WireDia.^2./4.*TLs.*CopperDensity;
    WeightSec_Insu = (pi.*(Sec_FullWireDia.^2 - Sec_WireDia.^2)./4).*TLs.*WireInsulationDensity;

    % Tape weight
    %---------------------------------------------
    
        
    if ~layerTapeUse
        Weight_InterlayerTape = zeros(size(Mlp));  % g
        V_tape = zeros(size(Mlp)); % g
        tapeMargin = 0;
    else
        tapeMargin = max(0.02*H, MinTapeMargin);
        % total build incl. tape (for average circumference)
        a1 = Mlp.*Pri_FullWireDia + max(Mlp-1,0).*numTapePerLayerPri.*kaptonThickness; % m
        a2 = Mls.*Sec_FullWireDia + max(Mls-1,0).*numTapePerLayerSec.*kaptonThickness; % m
    
        % base per-turn length near inner radius
        Lbase_p = zeros(size(Mlp));
        Lbase_s = zeros(size(Mls));
    
        if Winding_Pattern==1  % center-leg
            Lbase_p(isEE|isU)  = 2.*(PriW(isEE|isU) + PriH(isEE|isU) + 4.*CoreInsulationThickness(isEE|isU));
            Lbase_s(isEE|isU)  = Lbase_p(isEE|isU);
            Lbase_p(isER|isUR) = 2.*pi.*(PriW(isER|isUR)./2 + CoreInsulationThickness(isER|isUR));
            Lbase_s(isER|isUR) = Lbase_p(isER|isUR);
        else                   % double-leg
            Lbase_p(isEE|isU)  = 2.*(PriW(isEE|isU) + PriH(isEE|isU) + 4.*CoreInsulationThickness(isEE|isU));
            Lbase_p(isER|isUR) = 2.*pi.*(PriW(isER|isUR)./2 + CoreInsulationThickness(isER|isUR));
            Lbase_s(isEE)      = 2.*(SecW(isEE) + SecH(isEE) + 4.*CoreInsulationThickness(isEE));
            Lbase_s(isER)      = 2.*pi.*(sqrt(2).*PriW(isER)./2 + CoreInsulationThickness(isER));
            Lbase_s(isU)       = 2.*(PriW(isU) + PriH(isU) + 4.*CoreInsulationThickness(isU));
            Lbase_s(isUR)      = 2.*pi.*(PriW(isUR)./2 + CoreInsulationThickness(isUR));
        end
    
        % average circumference of tape wraps (≈ mid-build)
        Lavg_il_p = Lbase_p + pi.*(a1./2);
        Lavg_il_s = Lbase_s + pi.*(a2./2);
    
        % number of inter-layer boundaries
        nBndPri = max(Mlp-1,0);
        nBndSec = max(Mls-1,0);
    
        % total tape length (all boundaries × plies × overlap factor) (1.05 is 5% extra tape)
        L_tape_total = 1.05.*( nBndPri.*numTapePerLayerPri.*Lavg_il_p + nBndSec.*numTapePerLayerSec.*Lavg_il_s ); % m
    
        % tape width and volume
        w_tape = H + 2*tapeMargin;                 % m
        V_tape = kaptonThickness .* w_tape .* L_tape_total;  % m^3
    
        % mass in grams (KaptonDensity in g/m^3)
        Weight_InterlayerTape = kaptonDensity .* V_tape;  % g
    end
   
    % Compute total weight
    %---------------------------------------------
    TotalWeight = Wcore + WeightPri_copper + WeightSec_copper + WeightPri_Insu + ...
        WeightSec_Insu + WeightCore_Insu+Weight_InterlayerTape;

    % Calculate leakage inductance (not verified or used in this code)
    %--------------------------------------------------------------------
    
    % In this section, I need to include 187-190 calculations

    %XCp = 1/(2.*pi.*fs.*Cp)

    % Magnetizing inductance
    % Lm = ui.*u0.*Ac.*Np.^2./Le;
    % XLm = 2.*pi.*fs.*Lm;
    % 
    % WireInsulThickness = (Pri_FullWireDia-Pri_WireDia)./2;

    % Lleak = (u0.*Np.^2.*Tls)./(H-2.*WireInsulThickness).*...
    %     (WireInsulThickness+(Pri_PerLayer.*Pri_FullWireDia+ ...
    %     Sec_PerLayer.*Sec_FullWireDia)./3);

    % Unit length turn-to-tun capacitance following equation 8
    % Ctt = 0; % look at source

    % Secondary winding self capacitance
    % CparaSelf = ((Ns./Np).^2).*(Sec_PerLayer.*(Sec_PerLayer+1) ...
    %   .*(2.*Sec_PerLayer+1))./(6.*Sec_PerLayer.^2).*(4.*...
    %   (Mls-1))./(Mls.^2).*Tls.*Ctt;

    % Filter good designs
    %---------------------------------------------

    % Filter by B
    B_index = find(Bm < BSAT*BSAT_discountX);
    [Bmin,BminIndex] = min(Bm);
        
    % Filter by Temperature and Power Loss
    P_loss_index = find(Pcopper + Pcore <= Po./etaXfmer - Po);
    Tafterloss_index = find(Tafterloss <= TmaxX);
    Tmin_index = find(Tafterloss >= TminX);
    [Tminimum,TminValIndex] = min(Tafterloss);
    [Pmin,PminValIndex] = min(Pcopper+Pcore);

    % Filter by weight
    TotalWeight_index = find(TotalWeight < MaxWeightX);
    [WMin,WminValIndex] = min(TotalWeight);

    % Filter by packing factor min and max
    if Winding_Pattern == 2
        CopperPacking  = (((pi.*(Pri_WireDia.^2))./4).*Np + ((pi.*(Sec_WireDia.^2))./4).*(Ns./2))./ (H.*W);
        OverallPacking = (((pi.*(Pri_FullWireDia.^2))./4).*Np + ((pi.*(Sec_FullWireDia.^2))./4).*(Ns./2))./(H.*W);
    else 
        CopperPacking  = (((pi.*(Pri_WireDia.^2))./4).*Np + ((pi.*(Sec_WireDia.^2))./4).*Ns)./ (H.*W);
        OverallPacking = (((pi.*(Pri_FullWireDia.^2))./4).*Np + ((pi.*(Sec_FullWireDia.^2))./4).*Ns)./(H.*W);
    end
    if layerTapeUse
        % might be incorrect due to approx. as H*t 
        if Winding_Pattern == 2
            CopperPacking  = (((pi.*(Pri_WireDia.^2))./4).*Np + ((pi.*(Sec_WireDia.^2))./4).*(Ns./2))./ ((H-2.*tapeMargin).*W);
            OverallPacking = (((pi.*(Pri_FullWireDia.^2)/4).*Np + (pi.*(Sec_FullWireDia.^2)/4).*(Ns./2)) ...
                + (H.* ( max(Mlp-1,0).*numTapePerLayerPri.*kaptonThickness ...
                + max(Mls-1,0).*numTapePerLayerSec.*kaptonThickness))) ./ ((H-2.*tapeMargin).*W);
        else 
            CopperPacking  = (((pi.*(Pri_WireDia.^2))./4).*Np + ((pi.*(Sec_WireDia.^2))./4).*Ns)./ ((H-2.*tapeMargin).*W);
            OverallPacking = (((pi.*(Pri_FullWireDia.^2)/4).*Np + (pi.*(Sec_FullWireDia.^2)/4).*(Ns)) ...
                + (H.* ( max(Mlp-1,0).*numTapePerLayerPri.*kaptonThickness ...
                + max(Mls-1,0).*numTapePerLayerSec.*kaptonThickness))) ./ ((H-2.*tapeMargin).*W);
        end
    end
    OverallPackingmin_index = find(OverallPacking >= minpackingfactorX);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactorX);
    [PackingMin,PackingMinValIndex] = max(OverallPacking);
    [PackingMax,PackingMaxValIndex] = min(OverallPacking);

    if isempty(P_loss_index)
        fprintf("Transformer Bottlenecked. Min Power loss out of all candidates: %.2f Index: %d",Pmin,PminValIndex);
        y = zeros(1,39);
        return
    end
    if isempty(Tafterloss_index)
        fprintf("Transformer Bottlenecked. Min T out of all candidates: %.2f Index: %d",Tminimum,TminValIndex);
        y = zeros(1,39);
        return
    end
    if isempty(B_index)
        fprintf("Transformer Bottlenecked. Min B out of all candidates: %.2f Index: %d",Bmin,BminIndex);
        y = zeros(1,39);
        return
    end
    if isempty(TotalWeight_index)
        fprintf("Transformer Bottlenecked. Min Weight out of all candidates: %.2f Index: %d",WMin,WminValIndex);
        y = zeros(1,39);
        return
    end
    if isempty(OverallPackingmin_index)
        fprintf("Transformer Bottlenecked. Min packing factor out of all candidates: %.2f Index: %d",PackingMin,PackingMinValIndex);
        y = zeros(1,39);
        return
    end
    if isempty(OverallPackingmax_index)
        fprintf("Transformer Bottlenecked. Max packing factor out of all candidates: %.2f Index: %d",PackingMax,PackingMaxValIndex);
        y = zeros(1,39);
        return
    end

    % Determines viable turn amount based on window space
    %----------------------------------------------------

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
            (Pri_PerLayer.*Pri_FullWireDia <= H - 2.*CoreInsulationThickness-2*tapeMargin);
    
        EEER_center_sec_height_ok = is_EE_or_ER_core & ...
            (Sec_PerLayer.*Sec_FullWireDia <= H - 2.*CoreInsulationThickness-2*tapeMargin);
    
        EEER_center_width_ok = is_EE_or_ER_core & ...
            (Mlp.*Pri_FullWireDia+tTapePri+tTapeSec+(Mls.*Sec_FullWireDia) ...
             <= W - 3.*CoreInsulationThickness);
    
        meetsWindowFitConstraints = meetsWindowFitConstraints | ...
            (EEER_center_pri_height_ok & EEER_center_sec_height_ok & EEER_center_width_ok);
    end
    
    % ===== EE/ER cores: double-leg pattern =====
    if is_double_leg_pattern
        EEER_double_pri_height_ok = is_EE_or_ER_core & ...
            (Pri_PerLayer.*Pri_FullWireDia <= H - 2.*CoreInsulationThickness-2*tapeMargin);
    
        EEER_double_sec_height_ok = is_EE_or_ER_core & ...
            (Ns_group3.*Sec_FullWireDia <= H - 2.*CoreInsulationThickness-2*tapeMargin);
    
        EEER_double_width_ok = is_EE_or_ER_core & ...
            (Mlp.*Pri_FullWireDia +tTapePri+tTapeSec+ Ns_group2.*Sec_FullWireDia <= W - 3.*CoreInsulationThickness);
    
        meetsWindowFitConstraints = meetsWindowFitConstraints | ...
            (EEER_double_pri_height_ok & EEER_double_sec_height_ok & EEER_double_width_ok);
    end
    
    % ===== U/UR cores (table gives one set; apply for either selected pattern) =====
    % Height: Pri_perlayer*PriFull + Ns1*SecFull <= H - 3*TcoreInsu
    UUR_height_ok = is_U_or_UR_core & ...
        (Pri_PerLayer.*Pri_FullWireDia + Ns_group1.*Sec_FullWireDia <= H - 3.*CoreInsulationThickness-2*tapeMargin);
    
    % Width #1: Mlp*PriFull + Ns2*SecFull <= W - 3*TcoreInsu
    UUR_width1_ok = is_U_or_UR_core & ...
        (Mlp.*Pri_FullWireDia+tTapePri+Ns_group2.*Sec_FullWireDia <= W - 3.*CoreInsulationThickness);
    
    % Width #2: Mls*SecFull + Ns2*SecFull <= W - 2*TcoreInsu
    UUR_width2_ok = is_U_or_UR_core & ...
        (Mls.*Sec_FullWireDia+tTapeSec+Ns_group2.*Sec_FullWireDia ...
        <= W - 2.*CoreInsulationThickness);
    
    meetsWindowFitConstraints = meetsWindowFitConstraints | ...
        (UUR_height_ok & UUR_width1_ok & UUR_width2_ok);
    
    % Also require that the chosen layer/row plan can realize Ns overall
    meetsWindowFitConstraints = meetsWindowFitConstraints & (SupposeNs >= Ns);
    
    % Export indices used downstream
    WindowFit_index = find(meetsWindowFitConstraints);

    if isempty(WindowFit_index)
        fprintf("Transformer Bottlenecked. Max window area");
        y = zeros(1,39);
        return
    end

    % Final Filter
    %---------------------------------------------------------------------------------

    Index_Meet_All = intersect(B_index,P_loss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tmin_index);
    Index_Meet_All = intersect(Index_Meet_All,TotalWeight_index);
    Index_Meet_All = intersect(Index_Meet_All,OverallPackingmin_index);
    Index_Meet_All = intersect(Index_Meet_All,OverallPackingmax_index);
    Index_Meet_All = intersect(Index_Meet_All,WindowFit_index);

    % Sort by total weight and keep only the lightest one
    % [~,SortIndex] = sort(TotalWeight(Index_Meet_All));
    
    % Build Results Table 
    %% -----------------------------------------------------------------------------
    [~, SortIndex] = sort(TotalWeight(Index_Meet_All));
    if(length(SortIndex) >= 1)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1));

        V_cu     = ((WeightPri_copper+WeightSec_copper))/CopperDensity;
        V_insu   = ((WeightCore_Insu+WeightPri_Insu+WeightSec_Insu))/CoreInsulationDensity;
        Volume_m3 = Ve + V_cu + V_insu+V_tape;

        Sec_WireDiaMM = Sec_WireDia.*1000;
        Sec_WireAWG = -39*log(Sec_WireDiaMM./0.127)./log(92)+36;
        Pri_WireDiaMM = Pri_WireDia.*1000;
        Pri_WireAWG = -39*log(Pri_WireDiaMM./0.127)./log(92)+36;

        Design(:, 1)  = Po(TotalWeightSortIndex);
        Design(:, 2)  = Vppeak(TotalWeightSortIndex);
        Design(:, 3)  = Vspeak(TotalWeightSortIndex);
        Design(:, 4)  = fs(TotalWeightSortIndex);
        Design(:, 5)  = matno_record(TotalWeightSortIndex);
        Design(:, 6)  = matfs(TotalWeightSortIndex);
        Design(:, 7)  = Np(TotalWeightSortIndex);
        Design(:, 8)  = Ns(TotalWeightSortIndex);
        Design(:, 9)  = Bm(TotalWeightSortIndex);
        Design(:,10)  = Pri_WireAWG(TotalWeightSortIndex);
        Design(:,11)  = Pri_FullWireDia(TotalWeightSortIndex);
        Design(:,12)  = Sec_WireAWG(TotalWeightSortIndex);
        Design(:,13)  = Sec_FullWireDia(TotalWeightSortIndex);
        Design(:,14)  = Ippeak(TotalWeightSortIndex) ./ ...
                        (pi * Pri_Nstrands(TotalWeightSortIndex) .* Pri_ds(TotalWeightSortIndex).^2 / 4);
        Design(:,15)  = Ispeak(TotalWeightSortIndex) ./ ...
                        (pi * Sec_Nstrands(TotalWeightSortIndex) .* Sec_ds(TotalWeightSortIndex).^2 / 4);
        Design(:,16)  = Pri_Nstrands(TotalWeightSortIndex);
        Design(:,17)  = Sec_Nstrands(TotalWeightSortIndex);
        Design(:,18)  = Pri_PerLayer(TotalWeightSortIndex);
        Design(:,19)  = Mlp(TotalWeightSortIndex);
        Design(:,20)  = Sec_PerLayer(TotalWeightSortIndex);
        Design(:,21)  = Mls(TotalWeightSortIndex);
        Design(:,22)  = Ns_group1(TotalWeightSortIndex);
        Design(:,23)  = Ns_group2(TotalWeightSortIndex);
        Design(:,24)  = Ns_group3(TotalWeightSortIndex);
        Design(:,25)  = Ns_group4(TotalWeightSortIndex);
        Design(:,26)  = CopperPacking(TotalWeightSortIndex);
        Design(:,27)  = OverallPacking(TotalWeightSortIndex);
        Design(:,28)  = Pcore(TotalWeightSortIndex);
        Design(:,29)  = Pcopper(TotalWeightSortIndex);
        Design(:,30)  = Wcore(TotalWeightSortIndex);
        Design(:,31)  = WeightPri_copper(TotalWeightSortIndex);
        Design(:,32)  = WeightPri_Insu(TotalWeightSortIndex);
        Design(:,33)  = WeightSec_copper(TotalWeightSortIndex);
        Design(:,34)  = WeightSec_Insu(TotalWeightSortIndex);
        Design(:,35)  = WeightCore_Insu(TotalWeightSortIndex);
        Design(:,36)  = TotalWeight(TotalWeightSortIndex);
        Design(:,37)  = Tafterloss(TotalWeightSortIndex);
        Design(:,38)  = XcoreIndex(TotalWeightSortIndex);
        Design(:,39)  = Volume_m3(TotalWeightSortIndex);
        y = Design;
    else
        y = zeros(1,39);
        disp('Requirements not met. Filtered indexes == 0');
    end
end
end

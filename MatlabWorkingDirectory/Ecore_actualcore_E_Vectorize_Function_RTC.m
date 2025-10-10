function y = Ecore_actualcore_E_Vectorize_Function_RTC(...
    Vin_range,G_range,Po_range,Winding_Pattern, ...
    raw,raw1,raw2,raw3,raw4,raw5,raw6)

% Tunable Parameters
%% -------------------------------------------------------------------------------------


% Inductor parameters
%-------------------------------------------

% Lowest allowed inductor efficiency
etaInductor = 0.95;
% Max allowable temperature (C)
Tmax = 100;
% Min allowable temperature (C)
Tmin = 25;
% Max allowable current density in the wire (A/m^2)
Jwmax = 500*100*100;
% Minimal litz diameter one can get (m)
MinLitzDia = 0.05024/1000; %YAWG44, 0.0316 is AWG48, %0.03983 is
% Dielectric strength of the insulation material (V/m) , discount
dielectricstrength_insulation = 0.5*200*1000*100; %TEFLON
% minimal air gap (m)
mingap = 1e-4;
maxgap = 1e-3;
numGapsTested = 10;
% Minimum primary windings
MinWinding = 1;
% Maximum turns
MaxWinding = 100;
% Incremental winding
IncreN = 1;
% Maximum layer of winding
MaxMl = 10;
% Incremental layers
IncreMl = 1;
% Minimal wire diameter (m)
MinWireDia = 0.079/1000; %AWG28, 0.35 mm is AWG29, 0.079 is AWG40
% Maximum allowable weight (g)
MaxWeight = 100000;

% Electrical constants
%-------------------------------------------

% g/m3, density of copper
CopperDensity = 8.96*1000*1000;
% g/m3, density of core insulation materials
CoreInsulationDensity = 2.2*1000*1000; %TEFLON
% g/m3, density of wire insulation materials
WireInsulationDensity = 2.2*1000*1000; %IEFLON

% Discount factors
%--------------------------------------------

% Bmax discount factor
BSAT_discount = 0.85;
% Actual core loss is always higher than the calculated
CoreLossMultiple = 1.5;
% Maximum packing factor (copper area compared with total window area)
maxpackingfactor = 0.99;
% Minimum packing factor
minpackingfactor = 0.01;
% Winding factor of litz wire , assuming only 80% of wire size is copper
LitzFactor = 0.8;

% The variables and constants in this section are predefined. Recommend
% users not change them unless necessary

% Electrical constants. Normally there is no need to change
% ohm*m, resistivity of copper at 100C
rou = 2.3*1e-8;
% HA/m2, permeability of freespace
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
LCoreMU   = cell2mat(raw5(2:m1,3));
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
        FreqFlag(i) = 1;

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
% Column order:
% Ve(mm3), Ae (mm2) , Le (mm) , CoreShapeIndex, WindingPatternIndex, Window W (mm),...
% Half Window H (mm), Core W (mm),
% Half Core H (mm) , Core Thick (mm) , Center Leg Diameter (mm) , Total Window W (mm)
LCoreIndex = cell2mat(raw(2:m1, 1)) ;
LcoreVe = cell2mat(raw(2:m1,3))/(1000^3); % in m
LcoreAe = cell2mat(raw(2:m1,4))/(1000^2);
LcoreLe = cell2mat(raw(2:m1,5))/1000;
% This script is only for E core for now
LcoreCoreShapeIndex = cell2mat(raw(2:m1,6));
LcorePriW = cell2mat(raw(2:m1,8))/1000;
LcorePriH = cell2mat(raw(2:m1,9))/1000;
LcoreWindowW = cell2mat(raw(2:m1,12))/1000;
LcoreWindowH = cell2mat(raw(2:m1,13))/1000;

% DESIGN SWEEP
%% ------------------------------------------------------------------------

% This section lists all possible designs and prepare for the parfor loop
% in the next section

% Limit to core materials based on frequency
CoreMatIndexSweep = find(FreqFlag);

% Vectorize the design space
[Po, Vin ,G, matno_record , CoreIndex ,Np, Mlp, airgap] = ndgrid(Po_range ,Vin_range,...
G_range , CoreMatIndexSweep , LCoreIndex , MinWinding:IncreN:MaxWinding,1:IncreMl:MaxMl, linspace(mingap,maxgap,numGapsTested));

% Flatten
Po = reshape(Po,[] ,1);
Vin = reshape(Vin, [] , 1);
G = reshape(G,[] ,1);
matno_record = reshape(matno_record ,[] , 1);
Np = reshape(Np,[] ,1) ;
Mlp = reshape(Mlp,[] ,1);
airgap = reshape(airgap, [], 1);
CoreIndex = reshape(CoreIndex , [] , 1);


% Map CoreSize to actual size
Vcore = LcoreVe(CoreIndex); % in cm
Ac = LcoreAe(CoreIndex) ;
W = LcoreWindowW(CoreIndex);
H = LcoreWindowH(CoreIndex);
Le = LcoreLe(CoreIndex);
Center_L = LcorePriW(CoreIndex);
Center_T = LcorePriH (CoreIndex);
CoreShape = LcoreCoreShapeIndex(CoreIndex);

% Map material
ui = LCoreMU(matno_record);
BSAT = LCoreBSAT(matno_record);

% Inductance following RTC operating principles
Vo = Vin.*G;
L = u0*Ac.*Np.^2./(airgap + Le./ui);
% compute RTC boost timing and frequency
Ctot = 160e-12; % Assume using two GS66502T in parallel , 40 pF, two C3D1P7060 in parallel , 20 pF, 
% with some room

for xx = 1:1:length(L)
    sigma=(Vo(xx)-2*Vin(xx))./(Vo(xx)-Vin(xx));
    theta=acos(1-sigma);
    tring(xx)=(pi-theta).*sqrt(L(xx).*Ctot);
    Ilpeak(xx)=-(Vo(xx)-Vin(xx)).* sqrt(L(xx).*Ctot)./L(xx).* sin(theta);
    tlrise(xx)=L(xx).*abs(Ilpeak(xx))./Vin(xx);
    a=L(xx).*Vo(xx)./(Vo(xx)-Vin(xx));
    b=-2.*Po(xx).*Vo(xx).*L(xx)./Vin(xx)./(Vo(xx)-Vin(xx));
    c=-2.*Po(xx).*tring(xx)-2.*Po(xx).*tlrise(xx)-L(xx).*Ilpeak(xx).*Ilpeak(xx);
    d=-2.*Ctot.*Vo(xx).*Po(xx);
    p=[a b c d];
    I=roots(p);
    Ipeak(xx)=I(1);
end
tring = reshape(tring ,[] , 1);
Ilpeak = reshape(Ilpeak ,[] ,1);
tlrise = reshape(tlrise ,[] ,1);
Ipeak = reshape(Ipeak , [] ,1);

tfall = L.*Ipeak./(Vo-Vin);
trise = L.*Ipeak./Vin;
thold = Ctot.*Vo./Ipeak;
T = trise+thold+tfall+tring+tlrise;
fs=1./T;

% Not sure what this gives...
% ILave = (Ipeak.*(trise + thold + tfall)/2 + Ilpeak.*(tring + tlrise)/2)./T;

% Eliminate some elements based on dimension rule and BSAT rule
KeepAirGap = intersect(find(airgap >= mingap) , find(airgap <= 0.2*Le));
if isempty(KeepAirGap)
    warning("No feasible airgap: requested inductance larger than ungapped L0 for all cores/turns.");
end
ue = ui./(1+ui.*airgap./Le);
Bm_dummy = u0.*Np.*Ipeak./Le.*ue; %T

% Create filters for switching freq. and Bmax
Keep_Bmindex = find(Bm_dummy < BSAT*BSAT_discount);
Keep_fsindex = intersect(find(fs >= 1000000) ,find(fs <= 3000000));
KeepIndex = intersect(intersect(KeepAirGap,Keep_Bmindex),Keep_fsindex);

% Apply filters
Po = Po(KeepIndex);
Vin = Vin(KeepIndex);
Vo = Vo(KeepIndex);
Vinsulation_max = Vo;
matno_record = matno_record(KeepIndex);
BSAT = BSAT(KeepIndex);
CoreIndex = CoreIndex(KeepIndex);
Center_L = Center_L(KeepIndex);
Center_T = Center_T(KeepIndex);
CoreShape = CoreShape(KeepIndex);
H = H(KeepIndex);
W = W(KeepIndex);
Ac = Ac(KeepIndex);
Vcore = Vcore(KeepIndex);
airgap = airgap (KeepIndex);
Np = Np(KeepIndex);
Mlp = Mlp(KeepIndex);
L = L(KeepIndex);
fs = fs(KeepIndex);
Ipeak = Ipeak(KeepIndex);
Ilpeak = Ilpeak(KeepIndex);
Bm = Bm_dummy(KeepIndex);

% Find core loss property that's none zero around the required frequency for each design group
FsnoNonzero = F_atPv_500(matno_record,:) > 0;
FsnoIndex = abs(fs - F_atPv_500(matno_record ,:))./fs <= 0.4;
matfsIndex = FsnoNonzero.*FsnoIndex;
matfs = F_atPv_500(matno_record,:).*matfsIndex;
K1 = K1_range(matno_record,:).*matfsIndex*1000; %convert from mW/cm3 to W/m3
alpha = alpha_range(matno_record,:).*matfsIndex;
beta = beta_range(matno_record,:).*matfsIndex;
[rowIdcs , ~] = find(matfs > 0);

% So far, each row of the above represent one DESIGN POINT (that has one
% set of electrical requirements, one core size, one core material, one Np, Mlp and Mls)
% Each row of matfs, Ki, alpha and beta also correspond to each DESIGN POINT
% However, they have more than one non-zero columns because each material
% may have more than one loss data points in their datasheets around the required frequency

% We need to expand the design point to incorporate different loss data
% points for one material.

% Find the indices of unique values in rowIdcs
[UniqueRowIdcs , ~] = unique ( rowIdcs , 'rows') ;
ColDuplicate = sum(matfs(UniqueRowIdcs, :) ~= 0, 2);

% Repeat by the number of loss data of each design point
Po = repelem(Po(UniqueRowIdcs) ,ColDuplicate);
fs = repelem(fs(UniqueRowIdcs) , ColDuplicate);
Vin = repelem(Vin(UniqueRowIdcs) ,ColDuplicate);
Vo = repelem(Vo(UniqueRowIdcs), ColDuplicate);
Vinsulation_max = repelem(Vinsulation_max(UniqueRowIdcs) , ColDuplicate);
matno_record = repelem(matno_record (UniqueRowIdcs) , ColDuplicate);
BSAT = repelem (BSAT(UniqueRowIdcs) , ColDuplicate);
CoreIndex = repelem(CoreIndex(UniqueRowIdcs) ,ColDuplicate);
Center_L = repelem(Center_L(UniqueRowIdcs) ,ColDuplicate);
Center_T = repelem (Center_T(UniqueRowIdcs) ,ColDuplicate);
CoreShape = repelem(CoreShape(UniqueRowIdcs),ColDuplicate);
H = repelem(H(UniqueRowIdcs) ,ColDuplicate);
W = repelem (W( UniqueRowIdcs) ,ColDuplicate);
Ac = repelem(Ac(UniqueRowIdcs), ColDuplicate);
Vcore = repelem (Vcore (UniqueRowIdcs) , ColDuplicate);
airgap = repelem(airgap(UniqueRowIdcs) ,ColDuplicate);
Np = repelem(Np(UniqueRowIdcs) ,ColDuplicate);
Mlp = repelem(Mlp(UniqueRowIdcs) ,ColDuplicate);
L = repelem(L(UniqueRowIdcs) ,ColDuplicate);
Ipeak = repelem (Ipeak (UniqueRowIdcs) , ColDuplicate);
Ilpeak = repelem(Ilpeak(UniqueRowIdcs) ,ColDuplicate);
Bm = repelem(Bm(UniqueRowIdcs) , ColDuplicate);

% Reformat loss data into one non-zero vector
matfs = nonzeros(reshape(matfs(UniqueRowIdcs ,:)',[], 1));
K1 = nonzeros(reshape(K1(UniqueRowIdcs,:)',[],1));
beta = nonzeros(reshape(beta(UniqueRowIdcs,:)',[],1));
alpha = nonzeros(reshape(alpha(UniqueRowIdcs,:)',[],1));

% For the remaining indexes, compute the more detailed parameters 
% like losses, size, subsequent weight
% -------------------------------------------------------------------

if (isempty(Po))
    y = zeros(1,35);
else

    % Window area (m)
    Wa = H.*W;
    % Core weight (g)
    Wcore = Vcore.*CoreDensity(matno_record);
    % Core insulation thickness needed
    CoreInsulationThickness = Vinsulation_max./dielectricstrength_insulation;

    % Determine wire type, size, and num of strands if Litz
    % -----------------------------------------------

    % RMS current
    Iprms=Ipeak./sqrt(2);
    % AC skin depth
    skindepth=1./sqrt(pi.*fs.*u0./rou);
    % Area required of wire m^2
    Areq_p=Iprms./Jwmax;
    % solid equivalent diameter
    dsolid=2.*sqrt(Areq_p./pi);
    % Solid vs. litz
    useSolid=dsolid<=skindepth;
    % Litz diameter (only 1*skindepth here)
    dstrand_litz=max(MinLitzDia,skindepth);
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
        Pri_WireDia(idLitz)=2.*sqrt((Pri_Nstrands(idLitz).*Astrand(idLitz))./(pi.*LitzFactor));
    end
    % Strand diameter
    Pri_ds=max(MinWireDia,dsolid);
    Pri_ds(idLitz)=dstrand_litz(idLitz);
    % Full wire size with insulation
    Pri_FullWireDia=Pri_WireDia+2.*(Vin./dielectricstrength_insulation);

    A_pri_cu=(pi.*(Pri_WireDia.^2))./4;
    A_pri_full=(pi.*(Pri_FullWireDia.^2))./4;

    CopperPacking=(A_pri_cu.*Np)./(H.*W);
    OverallPacking=(A_pri_full.*Np)./(H.*W);

    % Computes mean length of turn for pri, accounting for geometry
    % and winding pattern
    %-----------------------------------------------------------------------------

    % Primary turns per layer
    Pri_PerLayer = floor(Np./Mlp);
    % Total length of windings (Assuming only E cores)
    TLp = Np.*2.*(Center_L + Center_T + 4*CoreInsulationThickness + 2.*Mlp.*Pri_FullWireDia);

    % Calculate Copper Loss & Core Loss
    %--------------------------------------------------------
    
    % Dowell Copper Loss
    PriKlayer = sqrt(pi.*Pri_Nstrands).*Pri_ds./2./(Pri_WireDia);
    Pri_xp = Pri_ds./2./skindepth.*sqrt(pi.*PriKlayer);
    Pri_Rdc = rou .* TLp ./ ( Pri_Nstrands .* (pi .* Pri_ds.^2 ./ 4) );
    Pri_Fr = Pri_xp.*((sinh(2.*Pri_xp) + sin(2.*Pri_xp))./(cosh(2.*Pri_xp) ...
        - cos(2.*Pri_xp)) + 2.*(Mlp.^2.*Pri_Nstrands - 1)./3.*(sinh(Pri_xp) ...
        - sin(Pri_xp))./(cosh(Pri_xp) + cos(Pri_xp)));
    Pri_Rac = Pri_Rdc.*Pri_Fr;

    % Standard Steinmetz core loss (W)
    Pcore = CoreLossMultiple.*Vcore.*K1.*fs.^alpha.*Bm.^beta;
    Pcopper = (Po./Vin).^2.*Pri_Rdc + (Ipeak - Ilpeak).^2/8.*Pri_Rac;
    
    % Calculate the temp rise
    %--------------------------------------------------------

    Rth = 16.31e-3.*(Ac.*Wa).^(-0.405);
    Tafterloss = Rth.*(Pcopper + Pcore) + 25;

    % Calculate the weight
    %----------------------------------------------------------------

    WeightPri_copper = pi.*Pri_WireDia.^2./4.*TLp.*CopperDensity;
    WeightPri_Insu = pi.*(Pri_FullWireDia.^2 - Pri_WireDia.^2)./4.*TLp.*WireInsulationDensity;
    WeightCore_Insu = (2.*H.*(Center_L + 2*Center_T) + 4.*W.*(Center_L + 2*Center_T) ...
        + H.*(2.*Center_L + 2*Center_T)).*CoreInsulationThickness.*CoreInsulationDensity;

    TotalWeight = Wcore + WeightPri_copper + WeightCore_Insu;

    % Filter the good designs
    %----------------------------------------------------------------

    % Bm is ~4.6e-4, very small
    % Core loss is ~0.1W
    B_index = find(Bm < BSAT*BSAT_discount);
    P_loss_index = find(Pcopper + Pcore <= Po*(1 - etaInductor));
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);
    TotalWeight_index = find(TotalWeight < MaxWeight);

    OverallPackingmin_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactor);
    Mlp_index = find(Mlp.*Pri_FullWireDia <= W - 2*CoreInsulationThickness);
    Pri_PerLayer_index = find(Pri_PerLayer.*Pri_FullWireDia < H - 2*CoreInsulationThickness);

    Index_Meet_All = intersect(B_index , P_loss_index);
    Index_Meet_All = intersect(Index_Meet_All, Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tmin_index);
    Index_Meet_All = intersect(Index_Meet_All, TotalWeight_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmin_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmax_index);
    Index_Meet_All = intersect(Index_Meet_All, Mlp_index);
    Index_Meet_All = intersect(Index_Meet_All, Pri_PerLayer_index);

    % Sort by total weight and keep only the lightest one
    %----------------------------------------------------------------
    [~ ,SortIndex] = sort (TotalWeight(Index_Meet_All));
    if (length(SortIndex) >= 1)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1:1));
        Volume = Vcore(TotalWeightSortIndex) ...
           + WeightPri_copper(TotalWeightSortIndex)./CopperDensity ...
           + WeightPri_Insu(TotalWeightSortIndex)./WireInsulationDensity ...
           + WeightCore_Insu(TotalWeightSortIndex)./CoreInsulationDensity;

        Design(:,1) = Po(TotalWeightSortIndex);
        Design(:,2) = Vin(TotalWeightSortIndex);
        Design(:,3) = Vo(TotalWeightSortIndex);
        Design(:,4) = Vinsulation_max(TotalWeightSortIndex);
        Design(:,5) = fs(TotalWeightSortIndex);
        Design(:,6) = matno_record (TotalWeightSortIndex);
        Design(:,7) = matfs(TotalWeightSortIndex);
        Design(:,8) = Center_L(TotalWeightSortIndex);
        Design(:,9) = Center_T(TotalWeightSortIndex);
        Design(:,10) = Ac(TotalWeightSortIndex);
        Design(:,11) = H(TotalWeightSortIndex);
        Design(:,12) = W(TotalWeightSortIndex);
        Design(:,13) = Np(TotalWeightSortIndex);
        Design(:,14) = Bm(TotalWeightSortIndex);
        Design(:,15) = Pri_WireDia(TotalWeightSortIndex);
        Design(:,16) = Pri_FullWireDia(TotalWeightSortIndex);
        Design(:,17) = Ipeak(TotalWeightSortIndex)./(pi*Pri_Nstrands ...
            (TotalWeightSortIndex).*Pri_ds(TotalWeightSortIndex).^2/4);
        Design(:,18) = Pri_Nstrands(TotalWeightSortIndex);
        Design(:,19) = Pri_PerLayer(TotalWeightSortIndex);
        Design(:,20) = Mlp(TotalWeightSortIndex);
        Design(:,21) = CopperPacking(TotalWeightSortIndex);
        Design(:,22) = OverallPacking(TotalWeightSortIndex);
        Design(:,23) = Pcore(TotalWeightSortIndex);
        Design(:,24) = Pcopper(TotalWeightSortIndex);
        Design(:,25) = Wcore(TotalWeightSortIndex);
        Design(:,26) = WeightPri_copper(TotalWeightSortIndex);
        Design(:,27) = WeightPri_Insu(TotalWeightSortIndex);
        Design(:,28) = WeightCore_Insu(TotalWeightSortIndex);
        Design(:,29) = TotalWeight(TotalWeightSortIndex);
        Design(:,30) = Tafterloss(TotalWeightSortIndex);
        Design(:,31) = L(TotalWeightSortIndex);
        Design(:,32) = airgap(TotalWeightSortIndex);
        Design(:,33) = CoreIndex(TotalWeightSortIndex);
        Design(:,34) = Volume;
        Design(:,35) = CoreShape(TotalWeightSortIndex);
        y = Design;
    else
        y = zeros(1,35);
        disp('Requirements not met. Filtered indexes == 0');
    end
end
end


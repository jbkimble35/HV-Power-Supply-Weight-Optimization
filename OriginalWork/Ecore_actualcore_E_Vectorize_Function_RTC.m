function y = Ecore_actualcore_E_Vectorize_Function_RTC(Date, RunNumber, Hypothesis, Notes,...
    Vin_range , G_range, Po_range , Vinsulation_max_range , Winding_Pattern)

CodeName = 'Ecore_actualcore_E_Vectorize_Function_RTC.m';

% Lowest allowed transformer efficiency
etaInductor = 0.95;
% Max allowable temperature (C)
Tmax = 90;
% Min allowable temperature (C)
Tmin = 25;
% Max allowable current density in the wire (A/m^2)
Jwmax = 500*100*100;
% Minimal litz diameter one can get (m)
MinLitzDia = 0.05024/1000; %YAWG44, 0.0316 is AWG48, %0.03983 is
% Dielectric strength of the insulation material (V/m) , discount
dielectricstrength_insulation = 0.5*200*1000*100; %IEFLON
% minimal air gap (m)
mingap = 100e-6;

MinWinding = 1;
% Maximum turns
MaxWinding = 50;
% Incremental winding
IncreN = 1;
% Maximum layer of winding
MaxMI = 10;
% Incremental layers
IncreMl = 1;
% Minimal wire diameter (m)
MinWireSize = 0.079/1000; %AWG28, 0.35 mm is AWG29, 0.079 is AWG40
% Maximum allowable weight (g)
MaxWeight = 1000;
% g/m3, density of the core
CoreDensity = 4.8*1000*1000;
% g/m3, density of copper
CopperDensity = 8.96*1000*1000;
% g/m3, density of core insulation materials
CoreInsulationDensity = 2.2*1000*1000; %TEFLON
% g/m3, density of wire insulation materials
WireInsulationDensity = 2.2*1000*1000; %IEFLON

% all discount factors
% Bmax discount factor
BSAT_discount = 0.75;
% Actual core loss is always higher than the calculated
CoreLossMultiple = 1;
% Maximum packing factor (copper area compared with total window area)
maxpackingfactor = 0.7;
% Minimum packing factor
minpackingfactor = 0.01;
% Winding factor of litz wire , assuming only 80% of wire size is copper
LitzFactor = 0.8;
% Weight of bobbin as a faction of the core insulation
BobbinWeightFactor = 0.5;

% The variables and constants in this section are predefined. Recommend
% users not change them unless necessary

% Simulation output variable section:
% Requirement file output
field1 = 'name';
value1_req = { 'Date' , 'RunNumber' , 'Hypothesis', 'Notes', 'CodeName',...
    'Vin_range(V)', 'G_range(V)', 'Po_range(W)',...
    'Vinsulation_max_range(V)', 'etaInductor', 'Tmin(C)' , 'Tmax(C)', 'Jwmax(A/m2)',...
    'MinLitzDia(m)' , 'dielectricstrength_insulation (V/m) ','min air gap (m)',...
    'Winding Pattern',...
    'MinWireSize(m)',...
    'MinWindingTurns', 'MaxWindingTurns', 'MaxNumberOfLayer', 'MaxWeight(g)',...
    'CoreDensity(g/m2) ', 'CopperDensity(g/m2) ' , 'CoreInsulationDensity (g/m2)',...
    'WireInsulationDensity(g/m2) ' , 'BSAT_discount', 'CoreLossMultiple',...
    'maxpackingfactor', 'minpackingfactor', 'LitzFactor' , 'BobbinWeightFactor'};

value1_design = {'Po(W) ', 'Vin(V)' , 'Va(V) ' , 'Vinsulation_ max (V) ' , 'fs (Hz) ' , 'matno',...
    'CoreMatFreq(Hz) ', 'CenterL (m) ','CenterT (m) ', 'CoreAc(m2) ','CoreWindowH(m) ',...
    'CoreWindowW(m) ','NumOfPri',...
    'BcoreDensity(T) ', 'WirePriDia(m) ' , 'WirePriFullDia(m)',...
    'WirePri_Idensity (A/m2) ',...
    'WirePriNstrands ' , 'WirePri_per_layer', 'WirePri_Nlayer' ,...
    'CopperPackingFactor ' , 'PackingFactor', 'LossCore(W)', ...
    'LossCopper(W)' , 'WeightCore(g) ' , 'WeightPri_copper (g) ' , 'WeightPri_Insu(g)',...
    'WeightCore_Insu(g)' , 'TotalWeight(g) ', 'TempAbsolute(C) ' , 'L' , 'airgap(m)',...
    'CoreIndex'};

field2 = 'data';
value2_req = {Date , RunNumber, Hypothesis , Notes ,CodeName, Vin_range , G_range , Po_range,...
    Vinsulation_max_range , etaInductor ,Tmin,Tmax,Jwmax,...
    MinLitzDia, dielectricstrength_insulation ,mingap...
    Winding_Pattern , MinWireSize ,...
    MinWinding , MaxWinding , MaxMI, MaxWeight, ...
    CoreDensity , CopperDensity , CoreInsulationDensity,...
    WireInsulationDensity , BSAT_discount , CoreLossMultiple,...
    maxpackingfactor , minpackingfactor , LitzFactor , BobbinWeightFactor};
value2_design = {[] ,[] ,[] ,[] ,[] ,[] , [] ,[] ,[] ,[] ,[] ,[] ,[] ,[],[],...
[] ,[] ,[] ,[] ,[] ,[] ,[] ,[] ,[], [],...
[] ,[] ,[], [] ,[], [] ,[], []};

global Requirement Design raw1 raw2 raw3 raw4 raw5 raw

% Save requirements
Requirement = struct(field1 , value1_req , field2 , value2_req);
% Save design results
Design = struct (field1, value1_design , field2 , value2_design);
% Electrical constants. Normally there is no need to change
% ohm*m, resistivity of copper at 100C
rou = 2.3*1e-8;
% /(ohm*m) , conductivity of copper
sigma = 1/rou;
% HA/m2, permeability of freespace
uO = 4*pi*10^(-7);
% F/m, permittivity of freespace
ebsl0 = 8.854*le-12;

%% MAIN BODY OF THE CODE STARTS FROM HERE
% Read the material properties to get the P at different B and F map from
% the CoreLossData.xlsx file
% This map will be used to calculate core loss later on
[m1,n1] = size(raw1);
XCoreMAT = raw1(2:m1,2);
XCoreFreq = cell2mat(raw1(2:m1,3:n1));
[m1,n1] = size(raw2);
XCoreBfield = cell2mat (raw2(2:m1,3:n1));
[m1,n1] = size(raw3);
XCorePloss = cell2mat(raw3(2:m1,3:n1));
[m1,n1] = size(raw4);
XCoreBSAT = cell2mat(raw4(2:m1,3:n1));
[m1,n1] = size(raw5) ;
XCoreMU = cell2mat (raw5(2:m1,3:n1));

% Constant
Pbar = 500; %500mW/cm3
PFfactor = 1;

% Draw out Pv plot vs B then interpolate
NoMat = m1-1;
Ball = 0.001:0.001:1;
MATcolorvector = rand(NoMat,3)
FreqFlag = zeros(size(1:1:NoMat));
for i = 1:1:NoMat
    DataSheetFreq = XCoreFreq(i, isnan(XCoreFreq(i,:)));
    NoFreq = length (DataSheetFreq)/2;
    colorvector = rand(NoFreq,3);
    for j = 1:1:NoFreq
        %Pv = ConstantA*Bfield+ConstantB
        ConstantA(i, j) = (log10(XCorePloss(i,2*j)) - log10(XCorePloss( i ,2*j -1)))/(log10(XCoreBfield(i,2*j)) - log10(XCoreBfield(i,2*j-1)));
        ConstantB(i, j) = log10(XCorePloss(i ,2*j)) - ConstantA(i,j)*log10(XCoreBfield(i,2*j));
        B_atPv_500(i, j) = 10^((log10(Pbar) - ConstantB(i,j))/ConstantA(i,j)); % in T
        F_atPv_500(i, j) = DataSheetFreq(2*j-1); % in Hz
        PF_atPv_500(i, j) = B_atPv_500(i, j)*F_atPv_500(i, j)&PFfactor;

        % if (abs(fs-range - F-atPv-500(i,j))./fs-range <= 0.4)
            FreqFlag(i) = 1;
        % end
        % Steinmetz
        if (j > 1)
            beta_range(i ,j) = log10(XCorePloss(i,2*j)/XCorePloss(i ,2*j-1))/log10(XCoreBfield(i ,2*j )/XCoreBfield(i ,2*j -1));
            %Third point
            XCorePloss_3rd(i, j) = 10.^(ConstantA(i, j -1)*log10(XCoreBfield(i,2*j)) + ConstantB(i ,j-1));
            alpha_range(i, j) = logl0(XCorePloss_3rd(i ,j)/XCorePloss(i,2*j))/log10(DataSheetFreq(2*j -3)/DataSheetFreq(2*j -1)); %(f2/fl)^alpha = P2/P1;
            K1_range(i,j) = XCorePloss(i ,2*j)/(XCoreBfield(i ,2*j)^ beta_range(i, j)) /(DataSheetFreq(2*j-1)^alpha_range(i, j)); %W/cm3
            %Repeat frequency 2's steinmetz parameter for frequency 1
            if (j == 2)
                beta_range(i ,j-1) = beta_range(i,j);
                alpha_range(i ,j-1)= alpha_range(i ,j);
                Kl_range(i ,j -1) = XCorePloss(i ,2*j-2)/(XCoreBfield(i ,2*j -2)^beta_range(i ,j -1))/(DataSheetFreq(2*j -3)^alpha_range(i ,j-1));
            end
        end
    end
end

% Core size
[m1,n1] = size(raw);
% Column order:
% Ve(mm3), Ae (mm2) , Le (mm) , CoreShapeIndex, WindingPatternIndex, Window W (nrn),...
% Half Window H (mm), Core W (mm),
% Half Core H (mm) , Core Thick (mm) , Center Leg Diameter (mm) , Total Window W (mm)
XCoreIndex = cell2mat(raw(3:m1, 1)) ;
XcoreVe = cell2mat(raw(3:m1,3))/(1000^3); % in m
XcoreAe = cell2mat(raw(3:m1,4))/(1000^2);
XcoreLe = cell2mat(raw(3:m1,5))/1000;
XcoreCoreShapeIndex = cell2mat(raw(3:m1,6));
XcorePriW = cell2mat(raw(3:m1,8))/1000;
XcorePriH = cell2mat(raw(3:m1,9))/1000;
XcoreSecW = cell2mat(raw(3:m1,10))/1000;
XcoreSecH = cell2mat(raw(3:m1,11))/1000;
XcoreWindowW = cell2mat(raw(3:m1,12))/1000;
XcoreWindowH = 2*cell2mat(raw(3:m1,13))/1000;

%
% DESIGN STARTS FROM HERE
% This section lists all possible designs and prepare for the parfor loop
% in the next section

% Limit to core materials based on frequency
CoreMatIndexSweep = find(FreqFlag);

% Vectorize the design space
[Po, Vin ,G, matno_record , CoreIndex ,Np, Mlp, airgap] = ndgrid(Po_range ,Vin_range,...
G_range , CoreMatIndexSweep , XCoreIndex , MinWinding:IncreN:MaxWinding,1:IncreMl:MaxMl, mingap:mingap:100*mingap);
Po = reshape(Po,[] ,1);
Vin = reshape(Vin , 1);
G = reshape(G,[] ,1);
matno_record = reshape(matno_record ,[] , 1);
Np = reshape(Np,[] ,1) ;
Mlp = reshape(Mlp,[] ,1);
airgap = reshape(airgap, [], 1);
CoreIndex = reshape(CoreIndex , [] , 1);

% Map CoreSize to actual size
Vcore = XcoreVe(CoreIndex); % in cm
Ac = XcoreAe(CoreIndex) ;
W = XcoreWindowW(Corelndex);
H = XcoreWindowH(CoreIndex);
Le = XcoreLe(CoreIndex);
Center_L = XcorePriW(CoreIndex);
Center_T = XcorePriH (CoreIndex);

% Map material
ui = XCoreMU(matno_record);
BSAT = XCoreBSAT(matno_record);

% Inductance following RTC operating principles
Vo = Vin.*G;
Vinsulation_max = Vo;
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
    c=-2.*Po(xx).*tring(xx)-2.*Po(xx).*tlrise(xx)-L(xx).*Ilpeak(xx)*Ilpeak(xx);
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
ILave = (Ipeak.*(trise + thold + tfall)/2 + Ilpeak.*(tring + tlrise)/2)./T;

% Eliminate some elements based on dimension rule and BSAT rule
KeepAirGap = intersect(find(airgap >= mingap) , find(airgap <= 0.2*Le));
ue = ui./(1+ui.*airgap./Le);
Bm_dummy = u0.*Np.*Ipeak./Le.*ue; %T
Keep_Bmindex = find(Bm_dummy < BSAT*BSAT_discount);
Keep_fsindex = intersect(find(fs >= 1000000) ,find(fs <= 3000000));
KeepIndex = intersect(intersect(KeepAirGap,Keep_Bmindex) ,Keep_fsindex);

Po = Po(KeepIndex);
Vin = Vin(KeepIndex);
G = G(KeepIndex);
Vo = Vo(KeepIndex);
Vinsulation_max = Vo;
matno_record = matno_record(KeepIndex);
ui = ui(KeepIndex);
BSAT = BSAT(KeepIndex);

CoreIndex = CoreIndex(KeepIndex);
Center_L = Center_L(KeepIndex);
Center_T = Center_T(KeepIndex);
H = H(KeepIndex);
W = W(KeepIndex);
Ac = Ac(KeepIndex);
Le = Le(Keeplndex);
Vcore = Vcore(KeepIndex);
airgap = airgap (KeepIndex);

Np = Np(KeepIndex);
Mlp = Mlp(Keeplndex);
L = L(KeepIndex);
fs = fs(KeepIndex);
Ipeak = Ipeak(KeepIndex);
Ilpeak = Ilpeak(KeepIndex);
trise = trise(KeepIndex);
Bm = Bm_dumrny(Keeplndex);

% Find core loss property that 's none zero around the required frequency for each design group
FsnoNonzero = F_atPv_500(matno_record,:) > 0;
FsnoIndex = abs(fs - F_atPv_500(matno_record ,:))./fs <= 0.4;
matfsIndex = FsnoNonzero.*FsnoIndex;
matfs = F_atPv_500(matno_record ,:).*matfsIndex;
K1 = Kl_range(matno_record ,:) .*matfsIndex*1000; %convert from mW/cm3 to W/m3
alpha = alpha_range(matno_record ,:) .*matfsIndex;
beta = beta_range(matno_record ,:) .* matfsIndex;
[rowIdcs , colIdcs] = find(matfs > 0);

% So far , each row of the above represent one DESIGN POINT (that has one
% set of electrical requirements , one core size , one core material , one Np, Mlp and MIs)
% Each row of matfs , Ki, alpha and beta also correspond to each DESIGN POINT
% However, they have more than one non-zero columns because each material
% may have more than one loss data points in their datasheets around the required frequency

% We need to expand the design point to incorporate different loss data
% points for one material.

% Find the indices of unique values in rowldcs
[UniqueRowIdcs , ind] = unique ( rowIdcs , 'rows') ;
ColDuplicate = sum(matfs(UniqueRowIdcs, :) == 0, 2);

% Repeat by the number of loss data of each design point
Po = repelem(Po(UniqueRowIdcs) ,ColDuplicate);
fs = repelem(fs(UniqueRowIdcs) , ColDuplicate);
Vin = repelem(Vin(UniqueRowIdcs) ,ColDuplicate);
Vo = repelem(Vo(UniqueRowIdcs), ColDuplicate);
Vinsulation_max = repelem(Vinsulation_max(UniqueRowIdcs) , ColDuplicate);
matno_record = repelem(matno_record (UniqueRowIdcs) , ColDuplicate);
ui = repelem(ui(UniqueRowIdcs) ,ColDuplicate);
BSAT = repelem (BSAT(UniqueRowIdcs) , ColDuplicate);

CoreIndex = repelem(CoreIndex(UniqueRowIdcs) ,ColDuplicate);
Center_L = repelem(Center_L(UniqueRowIdcs) ,ColDuplicate);
Center_T = repelem (Center_T(UniqueRowIdcs) ,ColDuplicate);
H = repelem(H(UniqueRowIdcs) ,ColDuplicate);
W = repelem (W( UniqueRowIdcs) ,ColDuplicate);
Ac = repelem(Ac(UniqueRowIdcs), ColDuplicate);
Le = repelem(Le(UniqueRowIdcs), ColDuplicate);
Vcore = repelem (Vcore (UniqueRowIdcs) , ColDuplicate);
airgap = repelem(airgap(UniqueRowIdcs) ,ColDuplicate);
L = repelem(L(UniqueRowIdcs) ,ColDuplicate);
Ipeak = repelem (Ipeak (UniqueRowIdcs) , ColDuplicate);
Ilpeak = repelem(Ilpeak(UniqueRowIdcs) ,ColDuplicate);
trise = repelem(trise(UniqueRowIdcs) ,ColDuplicate);
Bm = repelem(Bm(UniqueRowIdcs) , ColDuplicate);
% Reformat loss data into one non-zero vector
matfs = nonzeros(reshape(matfs(UniqueRowIdcs ,:)',[], 1));
K1 = nonzeros(reshape(K1(UniqueRowIdcs ,: )' , [] , 1));
beta = nonzeros(reshape(beta(UniqueRowIdcs ,:)' , [], 1));
alpha = nonzeros(reshape(alpha(UniqueRowIdcs ,:)' ,[],1));

%size (Po)
if (isempty(Po))
    y = 0;
else
    %Repeat elements by Primary Wire Number of Strands
    skindepth =1./sqrt(pi*fs*uO/rou);
    ds = max(skindepth , MinLitzDia*ones(size(skindepth))); % take the skin depth litz
    Pri_Nstrands = floor((Po*2/etaInductor./Vin/Jwmax)./(pi*ds.^2/4)) + 1;
    
    % Window area (m)
    Wa = H.*W;
    % Core weight (g)
    Wcore = Vcore .* CoreDensity;

    % Winding
    % Primary wire diameter (m)
    Pri_WireSize = sqrt(Pri_Nstrands.*pi .*ds.^2./4./LitzFactor./pi).*2;
    % Primary wire diameter (m) including the insulation layer
    Pri_FullWireSize = Pri_WireSize + (Vin./dielectricstrength_insulation).*2;

    CopperPacking = (pi.*Pri_WireSize.^2.*Np./4)./(H.*W);
    OverallPacking = (pi.*Pri_FullWireSize.^2.*Np./4)./(H.*W);

    % Winding structures
    % Core insulation thickness needed
    CoreInsulationThickness = Vinsulation_max./dielectricstrength_insulation;
    % Primary turns per layer
    Pri_PerLayer = floor(Np./Mlp);
    % Total length of windings
    TLp = Np.*2.*(Center_L + Center_T + 4*CorelnsulationThickness + 2.*Mlp.*Pri_FullWireSize);

    % Copper Loss parameters (Dowell)
    PriKlayer = sqrt(pi.*Pri_Nstrands).*ds./2./(Pri_WireSize);
    Pri_xp = ds./2./skindepth.*sqrt(pi.*PriKlayer);
    Pri_Rdc = rou.*TLp./(pi.*Pri_WireSize.^2./4);
    Pri_Fr = Pri_xp.*((sinh(2.*Pri_xp) + sin(2.*Pri_xp))./(cosh(2.*Pri_xp) - cos (2.*Pri_xp)) + 2.*(Mlp.^2.*Pri_Nstrands - 1)./3.*(sinh(Pri_xp) - sin(Pri_xp))./(cosh(Pri_xp) + cos(Pri_xp)));
    Pri_Rac = Pri_Rdc.*Pri_Fr;

    % Core loss and copper loss
    % Standard Steinmetz core loss (W)
    Pcore = CoreLossMultiple.*Vcore.*K1.*fs.^alpha.*Bm.^beta;
    Pcopper = (Po./Vin).^2.*Pri_Rdc + (Ipeak - Ilpeak).^2/8.*Pri_Rac;
    
    % Calculate the temp rise
    Rth = 16.31.*1e-3.*(Ac.*Wa).^(-0.405);
    Tafterloss = Rth.*(Pcopper + Pcore) + 25;

    % Calcualte the weight
    WeightPri_copper = pi.*Pri_WireSize.^2./4.*TLp.*CopperDensity;
    WeightPri_Insu = pi.*(Pri_FullWireSize.^2 - Pri_WireSize.^2)./4.*TLp.*WireInsulationDensity;
    WeightCore_Insu = (2.*H.*(Center_L + 2*Center_T) + 4.*W.*(Center_L + 2*Center_T) + H.*(2.*Center_L + 2*Center_T)).*CoreInsulationThickness.*CoreInsulationDensity;

    TotalWeight = Wcore + WeightPri_copper + WeightCore_Insu;

    % Filter the good designs
    B_index = find(Bm < BSAT*BSAT_discount);
    P_loss_index = find(Pcopper + Pcore <= Po*(1 - etaInductor));
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);
    TotalWeight_index = find(TotalWeight < MaxWeight);

    OverallPackingmin_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactor);
    Mlp_index = find(Mlp.*Pri_FullWireSize <= W - 2*CorelnsulationThickness);
    Pri_PerLayer_index = find(Pri_PerLayer.*PriFullWireSize < H - 2*CoreInsulationThickness);

    Index_Meet_All = intersect(B_index , P_loss_index);
    Index_Meet_All = intersect(Index_Meet_All, Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tmin_index);
    Index_Meet_All = intersect(Index_Meet_All, TotalWeight_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmin_index);
    Index_Meet_All = intersect(Index_Meet_All, OverallPackingmax_index);
    Index_Meet_All = intersect(Index_Meet_All, Mlp_index);
    Index_Meet_All = intersect(Index_Meet_All, Pri_PerLayer_index);

    % Sort by total weight and keep only the lightest five
    [WeightSort ,SortIndex] = sort (TotalWeight(Index_Meet_All));
    if (length(SortIndex) >= 5)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1:5));
        
        Design(1).data = Po(TotalWeightSortIndex);
        Design(2).data = Vin(TotalWeightSortIndex);
        Design(3).data = Vo(TotalWeightSortIndex) ;
        Design(4).data = Vinsulation_max (TotalWeightSortIndex);
        Design(5).data = fs(TotalWeightSortIndex) ;
        Design(6).data = matno_record (TotalWeightSortIndex);
        Design(7).data = matfs(TotalWeightSortIndex);
        Design(8).data = Center_L(TotalWeightSortIndex);
        Design(10).data = Ac(TotalWeightSortIndex);
        Design(11).data = H(TotalWeightSortIndex);
        Design(12).data = W(TotalWeightSortIndex);
        Design(13).data = Np(TotalWeightSortIndex);
        Design(14).data = Bm(TotalWeightSortIndex);
        Design(15).data = Pri_WireSize(TotalWeightSortIndex);
        Design(16).data = Pri_FullWireSize(TotalWeightSortIndex);
        Design(17).data = Ipeak(TotalWeightSortIndex)./(pi*Pri_Nstrands(TotalWeightSortIndex).*ds(TotalWeightSortIndex).^2/4);
        Design(18).data = Pri_Nstrands(TotalWeightSortIndex);
        Design(19).data = Pri_PerLayer(TotalWeightSortIndex);
        Design(20).data = Mlp(TotalWeightSortIndex);
        Design(21).data = CopperPacking(TotalWeightSortIndex);
        Design(22).data = OverallPacking(TotalWeightSortIndex);
        Design(23).data = Pcore(TotalWeightSortIndex);
        Design(24).data = Pcopper(TotalWeightSortIndex);
        Design(25).data = Wcore(TotalWeightSortIndex);
        Design(26).data = WeightPri_copper(TotalWeightSortIndex);
        Design(27).data = WeightPri_Insu(TotalWeightSortIndex) ;
        Design(28).data = WeightCore_Insu(TotalWeightSortIndex);
        Design(29).data = TotalWeight(TotalWeightSortIndex);
        Design(30).data = Tafterloss(TotalWeightSortIndex);
        Design(31).data = L(TotalWeightSortIndex);
        Design(32).data = airgap (TotalWeightSortIndex);
        Design(33).data = CoreIndex (TotalWeightSortIndex);
        y = 1;
    else
        y = 0;
    end
end
end


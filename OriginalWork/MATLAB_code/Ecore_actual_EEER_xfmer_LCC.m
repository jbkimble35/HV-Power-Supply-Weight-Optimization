function y = Ecore_actual_EEER_xfmer_LCC(raw,raw1,raw2,raw3,raw4,raw5,...
    Vppeak_range, Vspeak_range, Po_range, fs_range, Vinsulation_max_range, Winding_Pattern)

% Lowest allowed transformer efficiency
etaXfmer = 0.95;
% Max allowable temperature (C)
Tmax = 90;
% Min allowable temperature (C)
Tmin = 25;
% Max allowable current density in the wire (A/m^2)
Jwmax = 500*100*100;
% Minimal litz diameter one can get (m)
MinLitzDia = 0.07874/1000;% AWG 40. 0.05024/1000； %AWG44
% Dielectric strength of the insulation material (V/m), discount 50%
dielectricstrength_insulation = 0.5*200*1000*100; %IEFLON
% Minimum allowable core cross section radius (m)

MinPriWinding = 1;
% Maximum primary winding turns
MaxPriWinding = 100;
% Incremental pri winding
IncreNp = 1;
% Maximum layer of primary winding
MaxMlp = 10;
% Incremental pri layers
IncreMlp= 1;
% Maximum layer of secondary winding
MaxMls= 10;
% Incremental sec layers
IncreMls = 1;
% Minimal secondary wire diameter (m)
MinSecWireSize = 0.4/1000; %0.079 is AWG40, 0.4 mm is 178-5790
% Maximum allowable weight (g)
MaxWeight = 2000;
% g/m3, density of the core
CoreDensity = 4.8*1000*1000;
% g/m3, density of copper
CopperDensity = 8.96*1000*1000;
% g/m3, density of core insulation materials
CoreInsulationDensity = 2.2*1000*1000;%IEFLON
% g/m3, density of wire insulation materials
WireInsulationDensity = 2.2*1000*1000;%TEFLON

% all discount factors
% Bmax discount factor
BSAT_discount = 0.75;
% Actual core loss is always higher than the calculated.
CoreLossMultiple = 1.0;
% Maximum packing factor (copper area compared with total window area)
maxpackingfactor = 0.7;
% Minimum packing factor
minpackingfactor = 0.01;
% Winding factor of litz wire, assuming only 80% of wire size is copper
LitzFactor = 0.8;
% Weight of bobbin as a faction of the core insulation
BobbinWeightFactor = 0.5;

% Save design results
Design = zeros(1,43);
% Electrical constants. Normally there is no need to change
% ohm*m, resistivity of copper at 100C
rou = 2.3*1e-8;
% /(ohm*m), conductivity of copper
sigma = 1/rou;
% HA/m2, permeability of freespace
u0 = 4*pi*10^(-7);
% F/m, permittivity of freespace
ebs10 = 8.854*1e-12;
%% MAIN BODY OF THE CODE STARTS FROM HERE
% -----------------------------------------------------------------------
% Read the material properties to get the P at different B and F map from
% the CoreLossData.xlsx file
% This map will be used to calculate core loss later on
[m1,n1] = size(raw1);
XCoreMAT = raw1(2:m1,2);
XCoreFreq = cell2mat(raw1(2:m1,3:n1));
[m1,n1] = size(raw2);
XCoreBfield = cell2mat(raw2(2:m1,3:n1));
[m1,n1] = size(raw3);
XCorePloss = cell2mat(raw3(2:m1,3:n1));
[m1,n1] = size(raw4);
XCoreBSAT = cell2mat(raw4(2:m1,3:n1));
[m1,n1] = size(raw5);
XCoreMU = cell2mat(raw5(2:m1,3:n1));

%Constant
Pbar = 500; %500nW/cm3
PFfactor = 1; %

% Draw out Pv plot vs B then interpolate
NoMat = m1-1;
Ball = 0.001:0.001:1;
MATcolorvector = rand(NoMat,3);
FreqFlag = zeros(size(1:1:NoMat));
for i= 1:1:NoMat
    DataSheetFreq = XCoreFreq(i,~isnan(XCoreFreq(i,:)));
    NoFreq = length(DataSheetFreq)/2;
    colorvector = rand(NoFreq,3);
    for j= 1:1:NoFreq
        %Pv = ConstantA*Bfield+ConstantB
        ConstantA(i,j) = (log10(XCorePloss(i,2*j)) - log10(XCorePloss(i,2*j-1)))/( ...
            log10(XCoreBfield(i,2*j)) - log10(XCoreBfield(i,2*j-1)));
        ConstantB(i,j) = log10(XCorePloss(i,2*j)) - ConstantA(i ,j)*log10(XCoreBfield( ...
            i,2*j));
        B_atPv_500(i,j) = 10^((log10(Pbar) - ConstantB(i,j))/ConstantA(i,j)); % in T
        F_atPv_500(i,j) = DataSheetFreq(2*j-1); % in Hz
        PF_atPv_500(i,j) = B_atPv_500(i,j)*F_atPv_500(i,j)^PFfactor;
        
        if (abs(fs_range - F_atPv_500(i,j))./fs_range <= 0.4)
            FreqFlag(i) = 1;
        end
        %Steinmetz
        if (j > 1)
            beta_range(i,j) = log10(XCorePloss(i,2*j)/XCorePloss(i,2*j-1))/log10( ...
                XCoreBfield(i,2*j)/XCoreBfield(i,2*j-1));
            %Third point
            XCorePloss_3rd(i,j) = 10.^(ConstantA(i,j-1)*log10(XCoreBfield(i,2*j))+ ...
                ConstantB(i,j-1));
            alpha_range(i,j) = log10(XCorePloss_3rd(i,j)/XCorePloss(i,2*j))/log10( ...
                DataSheetFreq(2*j-3)/DataSheetFreq(2*j-l)); %(f2/f1)^alpha = P2/P1;
            K1_range(i,j) = XCorePloss(i,2*j)/(XCoreBfield(i,2*j)^beta_range(i,j))/( ...
                DataSheetFreq(2*j-1)^alpha_range(i,j)); %nW/cm3
            %Repeat frequency 2's steinmetz parameter for frequency 1
            if(j==2)
                beta_range(i,j-1) = beta_range(i,j);
                alpha_range(i,j-1) = alpha_range(i,j);
                K1_range(i,j-1) = XCorePloss(i,2*j-2)/(XCoreBfield(i,2*j-2)^ ...
                    beta_range(i,j-1))/(DataSheetFreq(2*j-3)^alpha_range(i,j-1));
            end
        end
    end
end
%Core size
[m1,n1] = size(raw);
TransformerCoreIndex = cell2mat(raw(3:m1,1));
XcoreVe = cell2mat(raw(3:m1,3))/(1000^3); % in m
XcoreAe = cell2mat(raw(3:m1,4))/(1000^2);
XcoreLe = cell2mat(raw(3:m1,5))/1000;
XcoreCoreShapeIndex = cell2mat(raw(3:m1,6));
XcorePriW = cell2mat(raw(3:m1,8))/1000;
XcorePriH = cell2mat(raw(3:m1,9)) /1000;
XcoreSecW = cell2mat(raw(3:m1,10)) /1000;
XcoreSecH = cell2mat(raw(3:m1,11)) /1000;
XcoreWindowW = cell2mat(raw(3:m1,12)) /1000;
XcoreWindowH = 2*cell2mat(raw(3:m1,13)) /1000;
ShuffleIndex = 1:1:length(TransformerCoreIndex);
%------------------------------------------------------------------------
% DESIGN STARTS FROM HERE
% Limit to core materials based on frequency
CoreMatIndexSweep = find(FreqFlag);
% Vectorize the design space

[Po, fs, Vppeak, Vspeak, Vinsulation_max,matno_record ,ShuffleXcoreIndex ,Np,Mlp,Mls] = ...
    ndgrid(Po_range , fs_range , Vppeak_range,...
    Vspeak_range, Vinsulation_max_range ,CoreMatIndexSweep, ShuffleIndex,...
    MinPriWinding:IncreNp:MaxPriWinding,1:IncreMlp:MaxMlp,1:IncreMls:MaxMls);

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
BSAT = XCoreBSAT( matno_record);
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

% Eliminate some elements based on dimension rule and BSAT rule
Bm_dummy = Vppeak/pi./fs./(2*Np.*Ac);
Keep_Bmindex = find(Bm_dummy < BSAT*BSAT_discount);
KeepIndex = Keep_Bmindex;

Po = Po(KeepIndex);
fs = fs(KeepIndex);
Vppeak = Vppeak(KeepIndex);
Vspeak = Vspeak(KeepIndex);
Vinsulation_max = Vinsulation_max(KeepIndex);
matno_record = matno_record(KeepIndex);
ui = ui(KeepIndex);
BSAT = BSAT( KeepIndex);
H = H(KeepIndex);
W = W( KeepIndex);
Ve = Ve(KeepIndex);
Ac = Ac(KeepIndex);
Le = Le(KeepIndex);
PriW = PriW(KeepIndex);
PriH = PriH(KeepIndex);
SecW = SecW(KeepIndex);
SecH = SecH(KeepIndex);
XcoreIndex = XcoreIndex(KeepIndex);
XcoreCoreShapeIndex = XcoreCoreShapeIndex(KeepIndex);

Np = Np(KeepIndex);
Mlp = Mlp(KeepIndex);
Mls = Mls(KeepIndex);

% Find core loss property that's none zero around the required frequency for each design group
FsnoNonzero = F_atPv_500(matno_record,:) > 0;
FsnoIndex = abs(fs - F_atPv_500(matno_record ,:))./fs <= 0.4;
matfsIndex = FsnoNonzero.*FsnoIndex;
matfs = F_atPv_500(matno_record,:).*matfsIndex;
K1 = K1_range(matno_record,:).*matfsIndex*1000; %convert from mW/cm3 to W/m3
alpha = alpha_range(matno_record,:).*matfsIndex;

beta = beta_range(matno_record,:).*matfsIndex;
[rowIdcs, colIdcs] = find(matfs > 0);
% So far, each row of the above represent one DESIGN POINT (that has one
% set of electrical requirements, one core size, one core material, one Np, Mlp and Mls);

% Each row of matfs, Kl, alpha and beta also correspond to each DESIGN POINT
% However, they have more than one non-zero columns because each material
% may have more than one loss data points in their datasheets around the required frequency

% We need to expand the design point to incorporate different loss data
% points for one material.

% Find the indices of unique values in rowIdcs
[UniqueRowIdcs, ind] = unique(rowIdcs,'rows');
ColDuplicate = sum(matfs(UniqueRowIdcs,:)~=0,2);
% Repeat by the number of loss data of each design point
Po = repelem(Po(UniqueRowIdcs),ColDuplicate);
fs = repelem(fs(UniqueRowIdcs),ColDuplicate);
Vppeak = repelem(Vppeak(UniqueRowIdcs) ,ColDuplicate);
Vspeak = repelem(Vspeak(UniqueRowIdcs),ColDuplicate);
Vinsulation_max = repelem(Vinsulation_max(UniqueRowIdcs) ,ColDuplicate);
matno_record = repelem(matno_record(UniqueRowIdcs) ,ColDuplicate);
ui = repelem(ui(UniqueRowIdcs), ColDuplicate);
BSAT = repelem (BSAT(UniqueRowIdcs),ColDuplicate);
H = repelem(H(UniqueRowIdcs),ColDuplicate);
W = repelem(W(UniqueRowIdcs) ,ColDuplicate);
Ve = repelem(Ve(UniqueRowIdcs), ColDuplicate);
Ac = repelem(Ac(UniqueRowIdcs),ColDuplicate);
Le = repelem(Le(UniqueRowIdcs),ColDuplicate);
XcoreIndex = repelem(XcoreIndex(UniqueRowldcs) ,ColDuplicate);
PriW = repelem(PriW(UniqueRowIdcs) ,ColDuplicate);
PriH = repelem(PriH(UniqueRowIdcs),ColDuplicate);
SecW = repelem(SecW(UniqueRowIdcs) ,ColDuplicate);
SecH = repelem(SecH(UniqueRowIdcs),ColDuplicate);
XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex(UniqueRowIdcs) ,ColDuplicate);

Np = repelem(Np(UniqueRowIdcs),ColDuplicate);
Mlp = repelem(Mlp(UniqueRowIdcs),ColDuplicate);
Mls = repelem(Mls(UniqueRowIdcs),ColDuplicate);
% Reformat loss data into one non-zero vector
matfs = nonzeros(reshape(matfs(UniqueRowldcs,:)' ,[],1));
K1 = nonzeros(reshape(K1(UniqueRowIdcs,:)' ,[],1));
beta = nonzeros(reshape(beta(UniqueRowIdcs,:)' ,[],1));
alpha = nonzeros(reshape(alpha(UniqueRowldcs ,:)' ,[],1));

if (isempty(Po))
    y = 0;
else
    %Repeat elements by Primary Wire Number of Strands
    skindepth = 1./sqrt(pi*fs*u0/rou);
    %ds = max(skindepth ,MinLitzDia*ones(size(skindepth))); % take the skin depth litz
    ds = MinLitzDia*ones(size(skindepth));
    MinPriNstrands = floor((Po*2/etaXfmer./Vppeak/Jwmax)./(pi*ds.^2/4))+ 1;
    MaxPriNstrands = floor((Po*2/etaXfmer./Vppeak/Jwmax*1.0)./(pi*ds.^2/4))+ 1;

    Po = repelem(Po,MaxPriNstrands - MinPriNstrands + 1);
    fs = repelem(fs ,MaxPriNstrands - MinPriNstrands + 1);
    Vppeak = repelem(Vppeak,MaxPriNstrands - MinPriNstrands + 1);
    Vspeak = repelem (Vspeak,MaxPriNstrands - MinPriNstrands + 1);
    Vinsulation_max = repelem(Vinsulation_max ,MaxPriNstrands - MinPriNstrands + 1);
    matno_record = repelem(matno_record,MaxPriNstrands - MinPriNstrands + 1);
    ui = repelem(ui,MaxPriNstrands - MinPriNstrands + 1);
    BSAT = repelem(BSAT,MaxPriNstrands - MinPriNstrands + 1);
    H = repelem(H, MaxPriNstrands - MinPriNstrands + 1);
    W = repelem(W,MaxPriNstrands - MinPriNstrands + 1);

    Ve = repelem(Ve,[MaxPriNstrands - MinPriNstrands + 1]);
    Ac = repelem(Ac,[MaxPriNstrands - MinPriNstrands + 1]);
    Le = repelem (Le, [MaxPriNstrands - MinPriNstrands + 1]);
    XcoreIndex = repelem(XcoreIndex,[MaxPriNstrands - MinPriNstrands + 1]);
    PriW = repelem(PriW,[MaxPriNstrands - MinPriNstrands + 1]);
    PriH = repelem(PriH,[MaxPriNstrands - MinPriNstrands + 1]);
    SecW = repelem (SecW, [MaxPriNstrands - MinPriNstrands + 1]);
    SecH = repelem(SecH,[MaxPriNstrands - MinPriNstrands + 1]);
    XcoreCoreShapeIndex = repelem(XcoreCoreShapeIndex,[MaxPriNstrands - MinPriNstrands + 1]);
    
    Np = repelem(Np,[MaxPriNstrands - MinPriNstrands + 1]);
    Mlp = repelem(Mlp,[MaxPriNstrands - MinPriNstrands + 1]);
    Mls = repelem(Mls,[ MaxPriNstrands - MinPriNstrands + 1]);
    matfs = repelem(matfs,[MaxPriNstrands - MinPriNstrands + 1]);
    K1 = repelem(K1,[MaxPriNstrands - MinPriNstrands + 1]);
    beta = repelem(beta,[MaxPriNstrands - MinPriNstrands + 1]);
    alpha = repelem(alpha ,[ MaxPriNstrands - MinPriNstrands + 1]);
    Pri_Nstrands = repmat((MinPriNstrands(1) :1: MaxPriNstrands(i))',length(MaxPriNstrands),1);

    %Repeat elements by Secondary Wire Number of Strands
    skindepth = 1./sqrt(pi*fs*u0/rou);
    ds = MinLitzDia*ones(size(skindepth));
    MinSecNstrands = 19*ones(size(skindepth));%178-5790 has 19 strands.
    MaxSecNstrands = 19*ones(size(skindepth));%178-5790 has 19 strands.

    Po = repelem(Po,[MaxSecNstrands - MinSecNstrands + 1]);
    fs = repelem(fs ,[MaxSecNstrands - MinSecNstrands + 1]);
    Vppeak = repelem(Vppeak,[MaxSecNstrands - MinSecNstrands + 1]);
    Vspeak = repelem(Vspeak ,[MaxSecNstrands - MinSecNstrands + 1]);
    Vinsulation_max = repelem(Vinsulation_max ,[MaxSecNstrands - MinSecNstrands + 1]);
    matno_record = repelem(matno_record ,[MaxSecNstrands - MinSecNstrands + 1]);
    ui = repelem(ui ,[MaxSecNstrands - MinSecNstrands + 1]);
    BSAT = repelem (BSAT,[MaxSecNstrands - MinSecNstrands + 1]);
    H = repelem (H,[ MaxSecNstrands - MinSecNstrands + 1]);
    W = repelem(W, [MaxSecNstrands - MinSecNstrands + 1]);
    Ve = repelem(Ve,[MaxSecNstrands - MinSecNstrands + 1]);
    Ac = repelem(Ac, [MaxSecNstrands - MinSecNstrands + 1]);
    Le = repelem (Le, [MaxSecNstrands - MinSecNstrands + 1]);
    XcoreIndex = repelem(XcoreIndex ,[MaxSecNstrands - MinSecNstrands + 1]);
    PriW = repelem(PriW,[MaxSecNstrands - MinSecNstrands + 1]);
    PriH = repelem(PriH,[MaxSecNstrands - MinSecNstrands + 1]);
    SecW = repelem (SecW, [MaxSecNstrands - MinSecNstrands + 1]);
    SecH = repelem (SecH, [MaxSecNstrands - MinSecNstrands + 1]);
    XcoreCoreShapeindex = repelem(XcoreCoreShapeIndex,[MaxSecNstrands - MinSecNstrands + 1]);
    Np = repelem(Np,[MaxSecNstrands - MinSecNstrands + 1]);
    Mlp = repelem(Mlp,[ MaxSecNstrands - MinSecNstrands + 1]);
    Mls = repelem (Mls ,[MaxSecNstrands - MinSecNstrands + 1]);
    matfs = repelem(matfs ,[MaxSecNstrands - MinSecNstrands + 1]);
    K1 = repelem (K1,[ MaxSecNstrands - MinSecNstrands + 1]);
    beta = repelem(beta ,[MaxSecNstrands - MinSecNstrands + 1]);
    alpha = repelem(alpha,[MaxSecNstrands - MinSecNstrands + 1]);
    Pri_Nstrands = repelem(Pri_Nstrands ,[MaxSecNstrands - MinSecNstrands + 1]);
    Sec_Nstrands = repmat((MinSecNstrands(1) :1: MaxSecNstrands(1))',length(MaxSecNstrands),1);
    
    % Recalculate several parameters
    k = Vspeak./ Vppeak;
    Ns = round(Np.*k)+1;
    % Primary current (A)
    Iprms = Po/etaXfmer./(Vppeak/sqrt(2));
    Ippeak = Iprms*sqrt(2);
    % Secondary current (A)
    Isrms = Po./(Vspeak/sqrt(2));
    Ispeak = Isrms*sqrt(2);

    skindepth = 1./sqrt(pi*fs*u0/rou);
    ds = MinLitzDia*ones(size(skindepth));
    % Calculate core loss (W)
    %-------------------------------------------------------------------
    % Peak flux, this corresponds to peak to peak flux density
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
    Wcore = Vcore.* CoreDensity;
    % Calculate Bmax (T)
    Bm = lamda./(2.*Np.*Ac);
    % Calculate core loss (W)
    Pcore = CorelossMultiple.Vcore.*Kl.fs.^alpha.*Bm.^beta;
    % Wire
    %----------------------------------------------------------------
    % Primary wire diameter (m)
    Pri_WireSize = sqrt(Pri_Nstrands.*pi.*ds.^2./4./LitzFactor./pi).*2;
    % Primary wire diameter (m) including the insulation layer
    Pri_FullWireSize = Pri_WireSize + (Vppeak./dielectricstrength_insulation).*2;
    % Secondary wire diameter (m)
    Sec_WireSize = sqrt(Sec_Nstrands.*pi.*ds.^2./4./LitzFactor./pi).*2;
    % Secondary wire diameter (m) including the insulation layer
    Sec_FullWireSize = Sec_WireSize + (Vspeak./ dielectricstrength_insulation /2).*2; 
        % ^for above line, 100% dielectric strength, similar with Rubadue data
    % For 178-5790 only
    Sec_FullWireSize = 1.016./1000*ones(size(Sec_WireSize));
    CopperPacking = (pi.*Pri_WireSize.^2.*Np./4 + pi.*Sec_WireSize.^2.*Ns/2./4)./(H.*W);
    OverallPacking = (pi.*Pri_FullWireSize.^2.*Np./4 + pi.*Sec_FullWireSize.^2.*Ns/2./4)./(H.*W);
    % Winding structures of each windings
    %----------------------------------------------------------------
    % Core insulation thickness needed
    CoreInsulationThickness = Vinsulation_max./dielectricstrength_insulation;
    % Total length of first layer of winding
    % Primary turns per layer
    Pri_PerLayer = floor(Np./Mlp);
    % Secondary turns per layer
    Sec_PerLayer = floor(Ns./Mls); %start with E core
    % Winding pattern of secondary
    % For ER, Cnly allow center leg winding
    switch Winding_Pattern
        case 1 %center leg
            % Secondary turns per layer
            Sec_PerLayer = floor(Ns./Mls);
            % Number of rows in section l, co-center with primary
            Ns_groupl = floor((H - 2*CoreInsulationThickness)./Sec_FullWireSize);
            Ns_group2 = zeros(size(Ns_group1));
            Ns_group3 = zeros(size (Ns_group1));
            Ns_group4 = zeros(size(Ns_group1));
            % Supposed secondary winding number if wind as mentioned above
            SupposeNs = Ns_groupl.*Mls;
            %% Total length of windings
            TLp = Np.*2.*(PriW + PriH + 4*CoreInsulationThickness + 2.*Mlp.* Pri_FullWireSize);
            TLs= Ns.*2.*(PriW + PriH + 4*Mlp.*Pri_FullWireSize + 8*CorelnsulationThickness + 2.*Mls.*Sec_FullWireSize);
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
            Sec_PerLayer = floor(Ns./2./Mis);
            % Winding on double leg
            Ns_group1 = zeros(size(Sec_PerLayer)); % does not wind on center leg
            Ns_group2 = floor((W - 3*CorelnsulationThickness - Mlp.* ...
            Pri_FullWireSize)./Sec_FullWireSize);
            Ns_group3 = floor((H - 2*CoreInsulationThickness - 2*Mls.* ...
            Sec_FullWireSize)./Sec_FullWireSize);
            Ns_group4 = Ns_group2;
            SupposeNs = Ns_group1.*Mls + Ns_group2.*Mls + Ns_group3.* ...
            Mls + Ns_group4.* Mls;
            %% Total length of windings
            TLp = Np.*2.*(PriW + PriH + 4*CoreInsulationThickness + 2.*Mlp.* ...
            Pri_FullWireSize);
            TLs = Ns.*2.*(SecW + SecH + 4*CoreInsulationThickness + 2.*Mls.* ...
            Sec_FullWireSize);
        otherwise
            disp('Wrong winding pattern');
    end
    % Calculate leakage inductance (not verified or used in this code)
    Lg = u0.*(W - Mls.*Sec_FullWireSize - Mlp.*Pri_FullWireSize).*SecH./H;
        %in Henry
    Xg = 2.*pi.*fs.*Lg;
    R_pri = Rload./Ns.^2;
    Lg_Lc_ratio = (W - 2.*CorelnsulationThickness - Mls.*Sec_FullWireSize - Mlp.* ...
    Pri_FullWireSize).*Le./ui./SecH./H;
    real_ratio = 1../(1 + Lg_Lc_ratio + Xg./ R_pri);
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
    Tafterloss = Rth.*(Pcopper + Pcore) + 25;
    % Calculate the weight
    %----------------------------------------------------------------
    WeightPri_copper = pi.*Pri_WireSize.^2./4.*TLp.* CopperDensity;
    WeightPri_Insu = pi.*(Pri_FullWireSize.^2 - Pri_WireSize.^2)./4.*TLp.* ...
    WireInsulationDensity;
    WeightSec_copper = pi.*Sec_WireSize.^2./4.*TLs.*CopperDensity;
    WeightSec_Insu = pi.*(Sec_FullWireSize.^2 - Sec_WireSize.^2)./4.*TLs.* ...
        WirelnsulationDensity;
    WeightCore_Insu =(2.*H.*(PriW + 2*PriH) + 4.*W.*(PriW + 2*PriH) + H.*(2*PriW + ...
        2*PriH)).*CoreInsulationThickness.*CoreInsulationDensity;
    % Recalculate XcoreCoreShapeIndex == 2, ER cores
    SelecIndex = find(XcoreCoreShapeIndex == 2);
    WeightCore_Insu(SelecIndex) - (sqrt(2)*pi*H(SelecIndex).*PriW(SelecIndex) + ...
        sqrt(2)*pi*2*W(SelecIndex).*PriW(SelecIndex) + H(SelecIndex)*pi.*PriW( ...
            SelecIndex)).*CoreInsulationThickness(SelecIndex).*CorelnsulationDensity;
    TotalWeight = Wcore + WeightPri_copper + WeightSec_copper + WeightPri_Insu + ...
        WeightSec_Insu + WeightCore_Insu;
    % Filter good designs
    %---------------------------------------------------------------
    B_index = find(Bm < BSAT*BSAT_discount);
    P_loss_index = find (Pcopper + Pcore <= Ploss_est);
    Tafterloss_index = find(Tafterloss <= Tmax);
    Tmin_index = find(Tafterloss >= Tmin);
    TotalWeight_index = find(TotalWeight < MaxWeight);

    OverallPackingmin_index = find(OverallPacking >= minpackingfactor);
    OverallPackingmax_index = find(OverallPacking <= maxpackingfactor);

    Mlp_index = find(Mlp.*Pri_FullWireSize <=W- 3*CoreInsulationThickness);
    Pri_PerLayer_index = find(Pri_PerLayer.*Pri_FullWireSize < H - 2* ...
    CorelnsulationThickness);
    % make sure pri and sec has enough insulation in between
    switch Winding_Pattern
        case 1 %center leg
            Mls_index = find(Mls.*Sec_FullWireSize + Mlp.*Pri_FullWireSize <=W- 3* ...
                CoreInsulationThickness);
        case 2 % double leg, Ns2 or Ns4 together with primary fits in the window 
            % width, Ns3 and primary also fits in the window width
            Mls_index = intersect(find(Ns_group2.*Sec_FullWireSize + Mlp.* ...
                Pri_FullWireSize <= W- 3*CoreInsulationThickness),...
                find(Ns_group3.*Sec_FullWireSize + Mlp.* Pri_FullWireSize <=W - 3* ...
                    CorelnsulationThickness));
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
    Index_Meet_All = intersect(Index_Meet.All,Tafterloss_index);
    Index_Meet_All = intersect(Index_Meet_All,Tmin_index);
    Index_Meet_All = intersect(Index_Meet.All,TotalWeight_index);
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
    [WeightSort,SortIndex] = sort(TotalWeight(Index_Meet_All));

    if(length(SortIndex) >= 1)
        TotalWeightSortIndex = Index_Meet_All(SortIndex(1:1));

        Design(:,1) = Po(TotalWeightSortIndex);
        Design(:,2) = Vppeak(TotalWeightSortIndex);
        Design(:,3) = Vspeak(TotalWeightSortIndex);
        Design(:,4) = Vinsulation_max(TotalWeightSortIndex);
        Design(:,5) = fs(TotalWeightSortIndex);
        Design(:,6) = matno_record(TotalWeightSortIndex);
        Design(:,7) = matfs(TotalWeightSortIndex);

        Design(:,8) = Ac(TotalWeightSortIndex);
        Design(:,9) = H(TotalWeightSortIndex);
        Design(:,10) =W(TotalWeightSortIndex);
        Design(:,11) = Np(TotalWeightSortIndex);
        Design(:,12) = Ns(TotalWeightSortIndex);
        Design(:,13) = real_ratio(TotalWeightSortIndex);
        Design(:,14) = Bm(TotalWeightSortIndex);

        Design(:,15) = Pri_WireSize(TotalWeightSortIndex);
        Design(:,16) = Pri_FullWireSize(TotalWeightSortIndex);
        Design(:,17) = Sec_WireSize(TotalWeightSortIndex);
        Design(:,18) = Sec_FullWireSize(TotalWeightSortIndex);
        Design(:,19) = Ippeak(TotalWeightSortIndex)./(pi+Pri_Nstrands( ...
            TotalWeightSortIndex).*ds(TotalWeightSortIndex).^2/4);
        Design(:,20) = Ispeak(TotalWeightSortIndex)./(pi*Sec_Nstrands( ...
            TotalWeightSortIndex).*ds(TotalWeightSortIndex).^2/4);

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
        y = Design;
    else
        y = zeros(1,43);
    end
end
end
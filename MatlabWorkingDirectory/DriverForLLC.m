clc, clf, clear

corelossfile = 'CoreLossData.xlsx';
coresizefile = 'CoreSizeData.xlsx';
coresizeSheetname = 'Ecore';

raw1 = readcell(corelossfile,'Sheet','Freq');
raw2 = readcell(corelossfile,'Sheet','Bfield');
raw3 = readcell(corelossfile,'Sheet','Ploss');
raw4 = readcell(corelossfile,'Sheet','BSAT');
raw5 = readcell(corelossfile,'Sheet','MU');
raw6 = readcell(corelossfile,'Sheet','Density');

% Ecore is the larger, perhaps inaccurate dataset, while ReviewedCores is a
% manually vetted selection of cores. OwnedCores is the 3 cores we own.
raw = readcell(coresizefile,'Sheet',coresizeSheetname);

%% Parameters to Adjust
%--------------------------------------------------------------------------

% Magnetizing inductance is assumed "large enough" as explained in the
% thesis, with A being the capacitance parallel vs. series ratio i.e.
% the inductance leakage vs. magnetizing ratio.

Date = '10-14-25';
% Quality factor
Q_range = 0.2:0.1:1;
% Resonant frequency
f0_range = 100000;
% frequency of the transformer
fs_range = f0_range;

% Capacitance ratio (inverse of inductance ratio) (shouldn't be lower than 0.1,
% since ZVS bandwidth becomes too small)
A_range = linspace(0.1,1,10);
% DC input voltage range (unipolar peak) (if Vppeak is the param. to select around,
% keep GT ~1, but optimal weight is usually achieved with tank gain of ~2)
Vin_range = 5;
% Peak of the output voltage that one hope to achieve (V)
% peak to peak is 2x this value
Vo_range = 20;
% Output power desired (W)
Po_range = 5;
% Turns ratio secondary/primary
K_range = 1:1:5;

% Insulation
%-------------------------------------------

% Use interlayer tape instead of full wire jacket ratings?
layerTapeUse = true;
enamelThickness = 20e-6;
kaptonDielStrength = 0.5*200e5; % V/m derated 50%
kaptonThickness = 60e-6; % m
MinTapeMargin = 5e-4;
kaptonDensity = 1.42e6; % g/m^3

% Core sheath insulation
CoreInsulationDensity = 2.2e6;       % g/m^3 (Teflon)
WireInsulationDensity = 2.2e6;       % g/m^3 (Teflon)
% Dielectric strength of core insulation material (V/m) 50% derated
dielectricstrength_insulation = 0.5 * 200e5;


% Inductor parameters
%-------------------------------------------

    % Lowest allowed inductor efficiency
    etaInductor = 0.98;
    % Max allowable temperature (C)
    TmaxL = 100;
    % Min allowable temperature (C)
    TminL = 25;
    % Maximum allowable weight (g)
    MaxWeightL = 1000;
    % Air gap (m)
    mingap = 0;
    
    % Winding and Wire Parameters
    %------------------------------------------
    
    % Minimum turns
    MinWindingL = 1;
    % Maximum turns
    MaxWindingL = 200;
    % Incremental winding
    IncreNL = 1;
    % Maximum layer of winding
    MaxMlL = 10;
    % Incremental layers. The layers of a transformer reference each wrap of
    % turns that fills the window height before moving on to the next level.
    % Once one layer fills, the next layer is wound on top, seperated by an
    % insulation layer.
    IncreMlL = 1;

    % Copper wire multiple to reduce resistive losses
    CuMultL = 1;
    
    % Discount factors
    %----------------------------------------
    
    % Bmax discount factor
    BSAT_discountL = 0.75;
    % Actual core loss is always higher than the calculated
    CoreLossMultipleL = 1;
    % Maximum packing factor (copper area compared with total window area)
    maxpackingfactorL = 0.7;
    % Minimum packing factor
    minpackingfactorL = 0.01;


% Transformer parameters
%-------------------------------------------
    
    % Minimum transformer efficiency
    etaXfmer = 0.95;
    % Max operating temp in Celsius
    TmaxX = 100;
    % Min operating temp in Celsius
    TminX = 25;
    % Minimum primary windings
    MinPriWindingX = 1;
    % Maximum primary windings
    MaxPriWindingX = 200;
    % Incremental primary winding
    IncreNpX = 1;
    % Maximum layer of primary winding
    MaxMlpX = 5;
    % Incremental layer of primary winding
    IncreMlpX = 1;
    % Maximum layer of secondary winding
    MaxMlsX = 25;
    % Incremental layer of secondary winding
    IncreMlsX = 1;
    % Max allowable transformer weight (g)
    MaxWeightX = 1000;

    % Deratings
    %------------------------------------------
    
    % Saturation flux density derating
    BSAT_discountX = 0.75;
    % Core loss multiplier
    CoreLossMultipleX = 1;
    maxpackingfactorX = 0.7;
    minpackingfactorX = 0.01;



% Winding factor of litz wire, assuming only 80% of wire size is copper
LitzFactor = 0.8;
% Minimal wire diameter (m)
MinWireDia = 0.25/1000; % AWG28, 0.35 mm is AWG29, 0.079 is AWG40
% Max allowable current density in the wire (A/m^2)
% 500A/cm^2 is the upper bound recommended, but without active cooling, and
% since the magnetics are thermally insulated, less is assumed
Jwmax = 3e6;
% Minimal litz diameter one can get (m)
MinLitzDia = 0.05024 / 1000; % AWG44 % 0.0316 is AWG48, 0.03983 is AWG46
% g/m^3, density of copper
CopperDensity = 8.96*1000*1000;
% Electrical constants. Normally there is no need to change
% ohm*m, resistivity of copper at 100C
rou = 2.3*1e-8;
% H/A·m^2, permeability of free space
u0 = 4*pi*10^(-7);


% Winding Pattern index: 1 indicates center leg winding, 2 indicates double
Winding_Pattern = 1;
% Hypothesis: record why you want to run the sim
Hypothesis ='';
% Notes: record any changes you made to the code
Notes ='';



%% File Output
%-------------------------------------------------------------------------------

% File output configuration
filename_xfmer = strcat(Date,'_','Xfmer.xlsx');
filename_inductor = strcat(Date,'_','Inductor.xlsx');
SheetNumber = 1;
Infosheetname = strcat('SimInfo', num2str(SheetNumber));
ResultDatasheetname = strcat('ResultsData', num2str(SheetNumber));

field1 = 'name';
value1_req = {'Date','Hypothesis','Notes',...
    'Q_range','fO_range','A_range','K_range',...
    'Vin_range','Vo_range','Po_range','fs_range','Winding_Pattern'};
field2 = 'data';
value2_req = {Date,Hypothesis, Notes,...
    Q_range, f0_range , A_range, K_range,...
    Vin_range, Vo_range, Po_range,fs_range , Winding_Pattern};
Requirement = struct(field1,value1_req,field2,value2_req);
fn   = fieldnames(Requirement);
vals = struct2cell(Requirement);

% Places input variable ranges and values in sheet named "SimInfo"
for i = 1:numel(vals)
    v = vals{i};
    if isnumeric(v) || islogical(v)
        if isscalar(v), vals{i} = v; else, vals{i} = mat2str(v); end
    elseif isstring(v) || ischar(v)
        vals{i} = char(v);
    else
        vals{i} = jsonencode(v);
    end
end
T = table(fn, vals, 'VariableNames', {'Field','Value'});
writetable(T, filename_xfmer, 'Sheet', Infosheetname, 'WriteVariableNames', true);

%% Calculations
%------------------------------------------------------------------------------

% First, resonant tank parameters are calculated

% Creates 4-D array of these 4 ranges. Each output is the size of all 4
% multiplied together. The reshape() just flattens each 4-D array into a
% column vector (turns ?x?x?x? into 4?x1)
[Q,f0,A,K] = ndgrid(Q_range, f0_range, A_range, K_range);
Q = reshape(Q,[],1);
f0 = reshape(f0,[],1);
A = reshape(A,[],1);
K = reshape(K,[],1);

% All of the following calculations are computed for each element in the
% matrix individually via the A.^B, A.*B, etc. operator. The sizes of A and B
% in A.^B must be equal or compatible. This allows for a large amount of
% independent values to be computed in a compact format.

% Equivalent resistance across secondary (from output p and output v)
Req = (Vo_range./sqrt(2)).^2./Po_range;
% Resonant tank reflected load resistance
RT = Req./K.^2;
% Series inductance of the resonant tank
Ls = RT./(2*pi.*f0.*Q);
% Series capacitance of LLC network
Cs = Q.*(A+1)./(A*2*pi.*f0.*RT);
% Parallel capacitance of LLC network
Cp = Q.*(A+1)./(2*pi.*f0.*RT);
% Transfer function gain factor; small-signal tank gain
GT = (4/pi)./(sqrt((1+A).^2.*(1-(fs_range./f0).^2).^2+1./Q.^2.*(fs_range./f0-A.*f0./((A+1).*fs_range)).^2));
% Maximum current through resonant tank
Imax = (Vin_range.*GT./RT).*sqrt(1+(fs_range./f0).^2.*Q.^2.*(A+1).^2);

% If effective gain GT.*K is within 20% of required gain, the design is
% acceptable. If not, index ignored. If all are ignored, error is thrown.
KeepIndex = intersect(find(GT.*K>=Vo_range./Vin_range),find(GT.*K<=1.2*Vo_range./Vin_range));
KeepIndex = intersect(KeepIndex,find(GT > 1));
if isempty(KeepIndex)
    error('Driver:NoCandidates', ...
          'No design points satisfy GT*K and GT>1. Adjust Vin/Vo/K/A/Q/f0.');
end

% Operating points not ignored are kept.
Q=Q(KeepIndex);
f0 = f0(KeepIndex);
A= A(KeepIndex);
K= K(KeepIndex);
RT = RT(KeepIndex);
Ls = Ls(KeepIndex);
Cs = Cs(KeepIndex);
Cp = Cp(KeepIndex);
GT = GT(KeepIndex);
Imax = Imax(KeepIndex);

% Loops over every row of the 4-D grid, with tic-toc measuring total
% runtime.
%-------------------------------------------
pcnt = 0.1;
tic
for i = 1:length(Q)
    
    % Peak voltage applied to primary from the input and resonant tank gain.
    Vpri = Vin_range.*GT(i);
    % Peak output voltage on the transformer after turn ratio K
    Vsec = Vin_range.*GT(i).*K(i);
    % Maximum insulation stress
    Vinsulation_max = Vsec;

    % Run Xfmer design, return design vector. All CoreLoss and CoreSize
    % data is passed, along with primary voltage, secondary voltage, output
    % power goal, switching frequency goal, max insulation stress, and^
    % winding pattern. None of the resonant tank sweep values are passed
    % here. Only the GT and K are relevant for the transformer design and
    % are what are being sweeped within this for loop through Vpri, Vinsulation_max,
    % and Vsec.
    SuceedX = Ecore_actual_EEER_xfmer_LCC(raw,raw1,raw2,raw3,raw4,raw5,raw6,...
        Vpri, Vsec, Po_range, fs_range, Vinsulation_max, Winding_Pattern,...
        layerTapeUse,enamelThickness,kaptonDielStrength,kaptonThickness,...
        MinTapeMargin,kaptonDensity,CoreInsulationDensity,WireInsulationDensity, ...
        dielectricstrength_insulation,etaXfmer,TmaxX,TminX,MinPriWindingX, ...
        MaxPriWindingX,IncreNpX,MaxMlpX,IncreMlpX,MaxMlsX,IncreMlsX,MaxWeightX, ...
        BSAT_discountX,CoreLossMultipleX,maxpackingfactorX,minpackingfactorX, ...
        LitzFactor,MinWireDia,Jwmax,MinLitzDia,CopperDensity,rou,u0);
    
    % Run Inductor design, return design table. All CoreLoss and CoreSize
    % data is passed, along with input voltage range (DC), output power
    % goal, switching frequency goal, winding pattern, and the resonant
    % tank sweep values.
    SuceedL = Ecore_actual_EEER_inductor_LCC(raw,raw1,raw2,raw3,raw4,raw5,raw6,...
        Vin_range,GT(i),Po_range,fs_range,Ls(i),Imax(i), Winding_Pattern,...
        Q(i), f0(i), A(i), K(i), RT(i), Ls(i), Cs(i), Cp(i), GT(i), ...
        layerTapeUse,enamelThickness,kaptonDielStrength,kaptonThickness,...
        MinTapeMargin,kaptonDensity,CoreInsulationDensity,WireInsulationDensity, ...
        dielectricstrength_insulation,etaInductor,TmaxL,TminL,MaxWeightL,mingap, ...
        MinWindingL,MaxWindingL,IncreNL,MaxMlL,IncreMlL,BSAT_discountL, ...
        CoreLossMultipleL,maxpackingfactorL,minpackingfactorL,CuMultL,...
        LitzFactor,MinWireDia,Jwmax,MinLitzDia,CopperDensity,rou,u0);
    
    % Successful result vector of many columns and 1 row for transformer and inductor are saved
    ResultX(i,:) = SuceedX;
    ResultL(i,:) = SuceedL;
    % Sliced variables in parallel loops allow this ResultX and ResultL to exist outside the
    % parfor loop.

    if i>=pcnt*length(Q)
        pcnt=pcnt+0.1;
        fprintf("%d Percent Complete \n",round(i*100/length(Q)));
    end
end
toc

% Results output
%-------------------------------------------

% Deletes rows of zeros, and then sorts by weight
XfmerDesignArray = ResultX(~all(ResultX == 0, 2), :);
XfmerDesignArray = sortrows(XfmerDesignArray,36,'ascend');

% Turns core geometry and material into their names from the sheet
freqTable = readcell(corelossfile,'Sheet','Freq');
sizeTable = readcell(coresizefile,'Sheet',coresizeSheetname);
matNames = freqTable(2:end,2);
geomNames = sizeTable(2:end,2);
geomIndexes = XfmerDesignArray(1:end,38);
matIndexes = XfmerDesignArray(1:end,5);
fullmatNames = matNames(matIndexes);
fullgeomNames = geomNames(geomIndexes);
XfmerDesignCellArr = num2cell(XfmerDesignArray);
XfmerDesignCellArr(:,5) = fullmatNames;
XfmerDesignCellArr(:,38) = fullgeomNames;

% Results for transformer and the column names are passed here.
XfmerDesignTable = cell2table(XfmerDesignCellArr,'VariableNames',{'Po_W','Vppeak_V',...
    'Vspeak_V','fs_Hz','Core Material','CoreMatFreq_Hz',...
    'NumOfPri','NumOfSec',...
    'BcoreDensity_T','WirePriDia_AWG','WirePriFullDia_m','WireSecDia_AWG',...
    'WireSecFullDia_m','WirePri_Idensity_A/m2','WireSecIdensity_A/m2',...
    'WirePriNstrands','WireSecNstrands','WirePri_per_layer','WirePri_Nlayer',...
    'WireSec_per_layer', 'WireSec_Nlayer','Ns1','Ns2','Ns3','Ns4','CopperPackingFactor',...
    'PackingFactor', 'LossCore_W','LossCopper_W' , 'WeightCore_g','WeightPri_copper_g',...
    'WeightPri_Insu_g', 'WeightSec_copper_g', 'WeightSec_Insu_g', 'WeightCore_Insu_g',...
    'TotalWeight_g', 'TempAbsolute_C','Core Geometry','Volume_m^3'});



% Results are written to excel file and sheet
xcelX = readcell(filename_xfmer,'Sheet',ResultDatasheetname);
[row,col] = size(xcelX);
writecell(repmat({''},row,col),filename_xfmer,'Sheet',ResultDatasheetname);
writetable(XfmerDesignTable,filename_xfmer,'Sheet',ResultDatasheetname);

if size(XfmerDesignTable,1)>=2
    weightX = XfmerDesignTable{2,36};
else
    weightX = 0;
end
fprintf("Transformer Weight is %.2f g",weightX);

% Deletes rows of zeros, and then sorts by weight
InductorDesignArray = ResultL(~all(ResultL == 0, 2), :);
InductorDesignArray = sortrows(InductorDesignArray,23,'ascend');

% Turns core geometry and material into their names from the sheet
freqTableL = readcell(corelossfile,'Sheet','Freq');
sizeTableL = readcell(coresizefile,'Sheet',coresizeSheetname);
matNamesL = freqTableL(2:end,2);
geomNamesL = sizeTableL(2:end,2);
geomIndexesL = InductorDesignArray(1:end,27);
matIndexesL = InductorDesignArray(1:end,5);
fullmatNamesL = matNamesL(matIndexesL);
fullgeomNamesL = geomNamesL(geomIndexesL);
InductorDesignCellArr = num2cell(InductorDesignArray);
InductorDesignCellArr(:,5) = fullmatNamesL;
InductorDesignCellArr(:,27) = fullgeomNamesL;

% Results for inductor and the column names are passed here.
InductorDesignTable = cell2table(InductorDesignCellArr,'VariableNames',{'PoW','Vin_V',...
    'Vpri_V','fs_Hz','Core Material','CoreMatFreq_Hz','NumOfPri','BcoreDensity_T', ...
    'WirePriDia_AWG','WirePriFullDia_m','WirePri_Idensity_Aperm2', ...
    'WirePriNstrands','WirePri_per_layer','WirePri_Nlayer',...
    'CopperPackingFactor', 'PackingFactor','LossCore_W',...
    'LossCopper_W','WeightCore_g', 'WeightPri_copper_g','WeightPri_Insu_g',...
    'WeightCore_Insu_g','TotalWeight_g','TempAbsolute_C','L', 'airgap_m', 'Core Geometry',...
    'Q','f0', 'A', 'K', 'RT', 'Ls', 'Cs', 'Cp', 'GT','Volume_m^3'});

xcelL = readcell(filename_inductor,'Sheet',ResultDatasheetname);
[row,col] = size(xcelL);
writecell(repmat({''},row,col),filename_inductor,'Sheet',ResultDatasheetname);

writetable(InductorDesignTable,filename_inductor,'Sheet',ResultDatasheetname);

if size(InductorDesignTable,1)>=2
    weightL = InductorDesignTable{2,23};
else 
    weightL = 0;
end
fprintf("Inductor Weight is %.2f g",weightL);

clc, clear

global Design_excel


% RTC test table is given in ComparisonTableRTC.xlsx

Date = 'xxx';
Vin_range = 170;
% Peak amplitude of the secondary voltage that one hope to achieve (V)
Vo_range = [362,363];
% Output power desired (W)
Po_range = 45;
% Winding Pattern Index: 1 inidcates center leg winding , 2 indicates double
Winding_Pattern = 1;
% Hypothesis: record why you want to run the sim
Hypothesis = 'xxx';
% Notes: record any changes you made to the code
Notes = 'xxx';

filename = strcat (Date, '_' , 'xxx.xlsx');
SheetNumber = 1;
Infosheetname = strcat('SimInfo', num2str(SheetNumber));
ResultDatasheetname = strcat( 'ResultsData' ,num2str(SheetNumber));
FBinvertersheetname = strcat('FBData' ,num2str(SheetNumber));
% RunNumber always starts with 1
RunNumber = 1;

% Read the core loss xlsx
corelossfile = 'CoreLossDataOLD.xlsx';
raw1 = readcell('CoreLossDataOLD.xlsx','Sheet','Freq');
raw2 = readcell('CoreLossDataOLD.xlsx','Sheet','Bfield');
raw3 = readcell('CoreLossDataOLD.xlsx','Sheet','Ploss');
raw4 = readcell('CoreLossDataOLD.xlsx','Sheet','BSAT');
raw5 = readcell('CoreLossDataOLD.xlsx','Sheet','MU');


% Read the core size xlsx
coresizefile = 'CoreSizeDataOLD.xlsx';
raw = readcell('CoreSizeDataOLD.xlsx','Sheet','Ecore');

% Read the FET and CAP xlsx
FETCapfile = 'MOSFETs and Capacitor Masses.xlsx';
FETsheetname = 'Component Masses_FETs';
CAPsheetname = 'Component Masses_Capacitors';

rawFET = readcell(FETCapfile,'Sheet',FETsheetname);
[m1,n1] = size(rawFET);
ColumnName = rawFET(1,:);

Vsw_term = 'voltage rating';  % V
Isw_term = 'current rating';  % A
Wsw_term = 'avg. mass (g)';    % g
Asw_term = 'Area (mm^2)';            % nm2

Wsw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Wsw_term))));
Isw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Isw_term))));
Vsw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Vsw_term))));
Asw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Asw_term))));

rawCAP = readcell(FETCapfile,'Sheet',CAPsheetname);
[m1,n1] = size(rawCAP);
ColumnName = rawCAP(1,:);

VCAP_term = 'Voltage (kV)';  % <V
CAP_term = 'Capacitance (uF)';  % uF
WCAP_term = 'Avg. mass (g)';  % g
ACAP_term = 'Area (mm^2)';  % nm2

WCAP = cell2mat(rawCAP(2:end, find(ismember(ColumnName, WCAP_term))));
CAP = cell2mat(rawCAP(2:end, find(ismember(ColumnName, CAP_term))));
VCAP = cell2mat(rawCAP(2:end, find(ismember(ColumnName, VCAP_term)))) * 1000;
ACAP = cell2mat(rawCAP(2:end, find(ismember(ColumnName, ACAP_term))));

PCBdensity = 0.0033; %g/mmmm
TotalCap = 1; %uF

tic
for i = 1:length(Vo_range)
    for Va = linspace(2*Vin_range,Vo_range(i),10)
        NumberOfStage = Vo_range(i)/Va;
        G = Va/Vin_range;
        Pa = Po_range/NumberOfStage;
        %YRun the Ecore Design
        Succeed = Ecore_actualcore_E_Vectorize_Function_RTC(Date,RunNumber, ...
            Hypothesis, Notes,Vin_range , G, Pa, Va, Winding_Pattern, ...
            raw,raw1,raw2,raw3,raw4,raw5);
        Result(i,:) = Succeed;
        NumberOfCopies(i) = NumberOfStage;
        
        % Run FB inverter selection
        VswIndex = find(Vsw >= 2*Va);
        SWparallel = floor(4.*Pa/Va./Isw(VswIndex))+1;
        SwitchWeight = 4*Wsw(VswIndex).*SWparallel;

        if size(SwitchWeight)==0
            warning('Switch not chosen');
        end

        [MinSwitchWeight(i), mi] = min(SwitchWeight);

        SWarea(i) = 4*Asw(VswIndex(mi))*SWparallel(mi);
        % Select Capacitors
        CAPseries = floor (2.*Va./VCAP) + 1;
        CAPparallel = floor(TotalCap./(CAP./CAPseries)) + 1;
        CAPWeight = WCAP.*CAPseries.*CAPparallel;
        [MinCapWeight(i),mi] = min(CAPWeight);
        CAParea(i) = ACAP(mi).*CAPseries(mi).*CAPparallel(mi)/2;
        % ^assume each cap space has 2 soldered on.
        % PCB weight
        PCBWeight(i) = 2*(SWarea(i) + CAParea(i))*PCBdensity;
        FBinverterWeight(i) = MinSwitchWeight(i) + ...
            MinCapWeight(i) + PCBWeight(i);
 
        FBdata = [NumberOfCopies(i),MinSwitchWeight(i),...
            MinCapWeight(i),PCBWeight(i),FBinverterWeight(i),...
            FBinverterWeight(i)*NumberOfCopies(i)];
        
        % Save Inductor Results

 
    end
end


names = {'Po(W)', 'Vin(V)' , 'Va(V)' , 'Vinsulation_ max (V)' , 'fs (Hz)' , 'matno',...
    'CoreMatFreq(Hz)', 'CenterL (m)','CenterT (m)', 'CoreAc(m2)','CoreWindowH(m)',...
    'CoreWindowW(m)','NumOfPri',...
    'BcoreDensity(T)', 'WirePriDia(m)' , 'WirePriFullDia(m)',...
    'WirePri_Idensity (A/m2)',...
    'WirePriNstrands' , 'WirePri_per_layer', 'WirePri_Nlayer' ,...
    'CopperPackingFactor' , 'PackingFactor', 'LossCore(W)', ...
    'LossCopper(W)' , 'WeightCore(g)' , 'WeightPri_copper (g)' , 'WeightPri_Insu(g)',...
    'WeightCore_Insu(g)' , 'TotalWeight(g)', 'TempAbsolute(C)' , 'L' , 'airgap(m)',...
    'CoreIndex'};
OutputTable = array2table(Design, 'VariableNames', names);
writetable(OutputTable, filename, 'Sheet', ResultDatasheetname);  % headers go to row 1


toc


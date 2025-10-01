clc, clear

% RTC test table is given in ComparisonTableRTC.xlsx

global FBdata Result

Date = '10-1-25';
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

PCBdensity = 0.0033; %g/mmmm
TotalCap = 1; %uF

filename = strcat (Date, '_' , 'RTC_inductor.xlsx');
SheetNumber = 1;
Infosheetname = strcat('SimInfo', num2str(SheetNumber));
ResultDatasheetname = strcat( 'ResultsData' ,num2str(SheetNumber));
FBinvertersheetname = strcat('FBData' ,num2str(SheetNumber));
% RunNumber always starts with 1
RunNumber = 1;

% Read the core loss xlsx
corelossfile = 'CoreLossData.xlsx';
raw1 = readcell('CoreLossData.xlsx','Sheet','Freq');
raw2 = readcell('CoreLossData.xlsx','Sheet','Bfield');
raw3 = readcell('CoreLossData.xlsx','Sheet','Ploss');
raw4 = readcell('CoreLossData.xlsx','Sheet','BSAT');
raw5 = readcell('CoreLossData.xlsx','Sheet','MU');
raw6 = readcell('CoreLossData.xlsx','Sheet','Density');

% Read the core size xlsx
coresizefile = 'CoreSizeData.xlsx';
raw = readcell('CoreSizeData.xlsx','Sheet','Ecore');

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

% Processes each output voltage target
%-----------------------------------------
pcnt = 0.1;
tic

Result = zeros(length(Vo_range),33);
FBdata = zeros(1,6);
for i = 1:length(Vo_range)

    % Per-stage output voltage chosen
    %---------------------------------

    for Va = linspace(2*Vin_range,Vo_range(i),10)
        
        % Per-stage spec.s

        NumberOfStage = Vo_range(i)/Va;
        G = Va/Vin_range;
        Pa = Po_range/NumberOfStage;

        % Run the Ecore Design for 1 stage

        Succeed = Ecore_actualcore_E_Vectorize_Function_RTC( ...
            Vin_range , G, Pa, Va, Winding_Pattern, ...
            raw,raw1,raw2,raw3,raw4,raw5,raw6);
        if Succeed(1,28)>Result(i,28)
            NumberOfCopies(i) = NumberOfStage;
            
            % Pick switches for full-bridge (per stage)
    
            VswIndex = find(Vsw >= 2*Va);
            SWparallel = floor(4.*Pa/Va./Isw(VswIndex))+1;
            SwitchWeight = 4*Wsw(VswIndex).*SWparallel;
            [MinSwitchWeight(i), mi] = min(SwitchWeight);
            SWarea(i) = 4*Asw(VswIndex(mi))*SWparallel(mi);
    
            % Select Capacitors (per stage)
    
            CAPseries = floor (2.*Va./VCAP) + 1;
            CAPparallel = floor(TotalCap./(CAP./CAPseries)) + 1;
            CAPWeight = WCAP.*CAPseries.*CAPparallel;
            [MinCapWeight(i),mi] = min(CAPWeight);
            CAParea(i) = ACAP(mi).*CAPseries(mi).*CAPparallel(mi)/2;
            % ^assume each cap space has 2 soldered on.
            % PCB weight
    
            % PCB weight (per stage) and total inverter weight
    
            PCBWeight(i) = 2*(SWarea(i) + CAParea(i))*PCBdensity;
            FBinverterWeight(i) = MinSwitchWeight(i) + ...
                MinCapWeight(i) + PCBWeight(i);
     
            FBdata(i,:) = [NumberOfCopies(i),MinSwitchWeight(i),...
                MinCapWeight(i),PCBWeight(i),FBinverterWeight(i),...
                FBinverterWeight(i)*NumberOfCopies(i)];

            Result(i,:)=Succeed;
        end
    end

    sortrows(FBdata,6,'ascend');
    OutputFB(i,:) = FBdata(1,:);

    sortrows(Result,28,"ascend");
    OutputL(i,:) = Result(1,:);

    if i>=pcnt*length(Vo_range)
        pcnt=pcnt+0.1;
        fprintf("%d Percent Complete \n",round(i*100/length(Vo_range)));
    end
end

OutputTableFB = array2table(FBdata,'VariableNames',{'Number of Stages', ...
    'Min Switch Weight (g)','Min Capacitor Weight (g)','PCB Weight (g)' ...
    ,'1 Stage Weight (g)','Total Weight (g)'});
arrFB = table2array(OutputTableFB);
OutputTableFB = OutputTableFB(~all(arrFB == 0, 2), :);
OutputTableFB = sortrows(OutputTableFB,'Total Weight (g)','ascend');

FBcelX = readcell(filename,'Sheet',FBinvertersheetname);
[row,col] = size(FBcelX);
writecell(repmat({''},row,col),filename_xfmer,'Sheet',FBinvertersheetname);
writetable(FBdata,filename,'Sheet',FBinvertersheetname);

OutputTableL = array2table(Result, 'VariableNames', {'Po(W)', 'Vin(V)' , ...
    'Va(V)' , 'Vinsulation_ max (V)' , 'fs (Hz)' , 'matno',...
    'CoreMatFreq(Hz)', 'CenterL (m)','CenterT (m)', 'CoreAc(m2)','CoreWindowH(m)',...
    'CoreWindowW(m)','NumOfPri','BcoreDensity(T)', 'WirePriDia(m)', ...
    'WirePriFullDia(m)','WirePri_Idensity (A/m2)','WirePriNstrands', ...
    'WirePri_per_layer', 'WirePri_Nlayer','CopperPackingFactor' , ...
    'PackingFactor', 'LossCore(W)','LossCopper(W)' , 'WeightCore(g)', ...
    'WeightPri_copper (g)' , 'WeightPri_Insu(g)','WeightCore_Insu(g)', ...
    'TotalWeight(g)', 'TempAbsolute(C)' , 'L' , 'airgap(m)','CoreIndex'});
arrL = table2array(OutputTableL);
OutputTableL = OutputTableL(~all(arrL == 0, 2), :);
OutputTableL = sortrows(OutputTableL,'TotalWeight(g)','ascend');

LcelX = readcell(filename,'Sheet',ResultDatasheetname);
[row,col] = size(LcelX);
writecell(repmat({''},row,col),filename_xfmer,'Sheet',ResultDatasheetname);
writetable(OutputTableL, filename, 'Sheet', ResultDatasheetname);

toc

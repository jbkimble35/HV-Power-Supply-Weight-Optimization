clc, clear all

global Requirement Design raw1 raw2 raw3 raw4 raw5 raw

Date = 'xxx';
Vin_range = 100;
% Peak amplitude of the secondary voltage that one hope to achieve (V)
Vo_range = 1000:1000:10000;
% Output power desired (W)
Po_range = 200;
% Winding Pattern Index: 1 inidcates center leg winding , 2 indicates double
Winding_Pattern = 1;
% Hypothesis: record why you want to run the sim
Hypothesis = ['xxx'];
% Notes: record any changes you made to the code
Notes = 'xxx';

filename = strcat (Date, '_' , 'xxx. xlsx ');
SheetNumber = 1;
Infosheetname = strcat('SimInfo', num2str(SheetNumber));
ResultDatasheetname = strcat( 'ResultsData' ,num2str(SheetNumber));
FBinvertersheetname = strcat('FBData' ,num2str(SheetNumber));
% RunNumber always starts with 1
RunNumber = 1;

% Read the core loss xlsx
corelossfile = 'CoreLossData.xlsx';
[num, txt ,raw1] = xlsread(corelossfile , 'Freq');
[num, txt ,raw2] = xlsread(corelossfile , 'Bfield');
[num, txt ,raw3] = xlsread(corelossfile , 'Ploss');
[num, txt ,raw4] = xlsread(corelossfile ,'BSAT');
[num, txt ,raw5] = xlsread(corelossfile ,'MU');

% Read the core size xlsx
coresizefile = 'CoreSizeData.xlsx';
[num, txt , raw] = xlsread(coresizefile , 'Ecore');

% Read the FET and CAP xlsx
FETCapfile = 'MOSFETs and Capacitor Masses.xlsx';
FETsheetname = 'Component Masses_FETs';
CAPsheetname = 'Component Masses_Capacitors';

%{
[num, txt ,rawFET] = xlsread(FETCapfile , FETsheetname);
[ml,n1] = size(rawFET);
ColumnName = rawFET(1,:);
Vsw_term = 'voltage rating '; %V
Isw_term = 'current rating '; %A
Wsw_term = 'avg. mass (g)'; %V
Asw_term = 'Area '; %nm2
Wsw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Wsw_term))));
Isw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Isw_term))));
Vsw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Vsw_term))));
Asw = cell2mat(rawFET(2:end, find(ismember(ColumnName, Asw_term))));

[num,txt ,rawCAP] = xlsread(FETCapfile ,CAPsheetname);
[ml,n1] = size(rawCAP);
ColumnName = rawCAP (1,:) ;
VCAP_term = 'Voltage (kV) '; %<V
CAP_term = 'Capacitance (uF) '; %uF
WCAP_term = 'avg. mass'; %g
ACAP_term = 'Area '; %nm2
WCAP = cell2mat(rawCAP(2:end, find(ismember(ColumnName, WCAP_term))));
CAP = cell2mat(rawCAP(2:end, find(ismember(ColumnName, CAP_term))));
VCAP = cell2mat(rawCAP(2:end, find(ismember (ColumnName, VCAP_term)))) *1000;
ACAP = cell2mat(rawCAP(2:end, find(ismember(ColumnName, Asw_term))));
%}
PCBdensity = 0.0033; %g/mmmm
TotalCap = 1; %uF

tic
for i = 1:1: length(Vo_range)
    for Va = linspace(2*Vin_range,Vo_range(i) ,10)
        NumberOfStage = Vo_range(i)/Va;
        G = Va/Vin_range;
        Pa = Po_range/NumberOfStage;
        
        %YRun the Ecore Design
        Suceed = Ecore_actualcore_E_Vectorize_Function_RTC(Date, RunNumber, Hypothesis, Notes,...
            Vin_range , G, Pa, Va, Winding_Pattern);
        if (Suceed == 1)
            NumberOfCopies(i) = NumberOfStage;
            
            %{
            % Run FB inverter selection
            VswIndex = find(Vsw >= 2*Va);
            SWparallel = floor(4.*Pa/Va./Isw(VswIndex))+1;
            SwitchWeight = 4*Wsw(Vswlndex).* SWparallel;
            [MinSwitchWeight(i), mi] = min(SwitchWeight);
            SWarea(i) = 4*Asw(VswIndex(mi))*SWparallel(mi);
            % Select Capacitors
            CAPseries = floor (2.*Va./VCAP) + 1;
            CAPparallel = floor(TotalCap./(CAP./CAPseries)) + 1;
            CAPWeight = WCAP.*CAPseries.*CAPparallel;
            [MinCapWeight(i),mi] = min(CAPWeight);
            CAParea(i) = ACAP(mi).*CAPseries(mi).*CAPparallel(mi)/2; %assume each cap space has 2 soldered on.
            % PCB weight
            PCBWeight(i) = 2*(SWarea(i) + CAParea(i))*PCBdensity;
            FBinverterWeight(i) = MinSwitchWeight(i) + MinCapWeight(i) + PCBWeight(i);
     
            FBdata = [NumberOfCopies(i),MinSwitchWeight(i),MinCapWeight(i),PCBWeight(i),FBinverterWeight(i),FBinverterWeight(i)*NumberOfCopies(i)];
            %}
            % Save Inductor Results
            Design_excel = squeeze(struct2cell(Design));
            [row,col] = size(Design_excel(2,:));
            Design_num = cell2mat (Design_excel(2,:));
            if (length(Design_excel (1,:) ) >= 26)
                CharName = strcat('A' ,char(64+length(Design_excel (1,:))-26));
                if (length(Design_excel(1,:)) >= 52)
                    disp ('Too Many Variables');
                end
            else
                CharName = char(64+length( Design_excel(1,:))-26);
            end
            xlswrite(filename , Design_num , ResultDatasheetname , strcat ('A', num2str(5*RunNumber -3),': ',CharName, num2str(5*RunNumber + 1)));
            % Save FB results
            xlswrite(filename ,[FBdata; FBdata; FBdata; FBdata; FBdata], FBinvertersheetname, strcat('A', num2str(5*RunNumber -3), ': ' , 'F', num2str(5*RunNumber + 1)));
            RunNumber = RunNumber + 1;
        end
    end
end
%%
% Save Requirements
Requirement_excel = squeeze(struct2cell(Requirement))';
xlswrite(filename ,Requirement_excel ,Infosheetname);
% Format Excel for easy to read
xlswrite(filename, Design_excel (1,:) ,ResultDatasheetname);
exl = actxserver( 'excel.application');

savefile = 'xxx';
exlWkbk = exl.Workbooks.Open(strcat(savefile, filename));
Worksheet = exlWkbk.Sheets.Item(3) ;
Worksheet.Rows.Item (1).RowHeight = 60;
cells = Worksheet.Range(strcat('A1: ',CharName, '1'));
cells.ColumnWidth = 6;
set(cells.Font, 'Bold', true);
exlWkbk.Save
exlWkbk.Close
exl.Quit
toc
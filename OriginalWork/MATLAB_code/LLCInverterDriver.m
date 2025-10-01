clc, clf, clear all
close all

% Read the core loss xlsx
corelossfile = 'CoreLossData.xlsx';
[num,txt,raw1] = xlsread(corelossfile,'Freq');
[num,txt,raw2] = xlsread(corelossfile,'Bfield');
[num,txt,raw3] = xlsread(corelossfile,'Ploss');
[num,txt,raw4] = xlsread(corelossfile,'BSAT');
[num,txt,raw5] = xlsread(corelossfile,'MU');

% Read the core size xlsx
coresizefile = 'CoreSizeData.xlsx';
[num,txt,raw] = xlsread(coresizefile,'Ecore');

Date = 'xxx';
% Quality factor
Q_range  = 0.1:0.1:10;
% Natural frequency
f0_range = 480000;
% Capacitance ratio
A_range  = 0.1:0.1:10;

% Turns ratio
K_range = 5:1:30;
% DC input voltage range
Vin_range = 200;
% Peak amplitude of the output voltage that one hope to achieve (V)
Vo_range = 7500;
% Output power desired (W)
Po_range = 750;
% frequency of the transformer
fs_range = 500000;

% Winding Pattern Index: 1 indicates center leg winding, 2 indicates double
Winding_Pattern = 1;
% Hypothesis: record why you want to run the sim
Hypothesis = [' '];
% Notes: record any changes you made to the code
Notes = '';

filename_xfmer    = strcat(Date,'_','xxx.xlsx');
filename_inductor = strcat(Date,'_','xxx.xlsx');
SheetNumber = 1;
Infosheetname = strcat('SimInfo',num2str(SheetNumber));
ResultDatasheetname = strcat('ResultsData',num2str(SheetNumber));

[Q,f0,A,K] = ndgrid(Q_range,f0_range,A_range,K_range);
Q  = reshape(Q ,[],1);
f0 = reshape(f0,[],1);
A  = reshape(A ,[],1);
K  = reshape(K ,[],1);

RT = 8/pi^2*40000^2/700/6/6./K.^2;
Ls = RT./(2*pi.*f0.*Q);
Cs = Q.*(A+1)./(A*2*pi.*f0.*RT);
Cp = Q.*(A+1)./(2*pi.*f0.*RT);
GT = 4/pi./(sqrt((1+A).^2.*(1-(fs_range./f0).^2).^2 + 1./Q.^2.*(fs_range./f0 - A.*f0./((A+1).*fs_range)).^2));
Imax = Vin_range.*GT./RT.*sqrt(1 + (fs_range./f0).^2.*Q.^2.*(A + 1).^2);

KeepIndex = intersect(find(GT.*K >= Vo_range/Vin_range), find(GT.*K <= 1.2*Vo_range/Vin_range));
KeepIndex = intersect(KeepIndex, find(GT > 1));

Q=Q(KeepIndex);
f0=f0(KeepIndex);
A=A(KeepIndex);
K=K(KeepIndex);
RT=RT(KeepIndex);
Ls=Ls(KeepIndex);
Cs=Cs(KeepIndex);
Cp=Cp(KeepIndex);
GT=GT(KeepIndex);
Imax=Imax(KeepIndex);

tic
parfor i = 1:1:length(Q)
    % Intermediate voltage
    Vpri = Vin_range.*GT(i);
    % Output voltage on the transformer
    Vsec = Vin_range.*GT(i).*K(i);
    Vinsulation_max = Vsec;

    % Run Xfmer design
    SucceedX = Ecore_actual_EEER_xfmer_LCC(raw,raw1,raw2,raw3,raw4,raw5,...
        Vpri, Vsec, Po_range, fs_range, Vinsulation_max, Winding_Pattern);

    SucceedL = Ecore_actual_EEER_inductor_LCC(raw,raw1,raw2,raw3,raw4,raw5,...
        Vin_range, GT(i), Po_range, fs_range, Ls(i), Imax(i), Winding_Pattern,...
        Q(i), f0(i), A(i), K(i), RT(i), Ls(i), Cs(i), Cp(i), GT(i));

    ResultX(i,:) = SucceedX;
    ResultL(i,:) = SucceedL;
end
toc

XfmerDesignTable = array2table(ResultX,'VariableNames',{'Po_W','Vppeak_V','Vspeak_V',...
    'Vinsulation_max_V','fs_Hz','matno','CoreMatFreq_Hz','CoreAc_m2','CoreWindowH_m',...
    'CoreWindowW_m','NumOfPri','NumOfSec','RealConversion','BcoreDensity_T','WirePriDia_m',...
    'WirePriFullDia_m','WireSecDia_m','WireSecFullDia_m','WirePri_Idensity_Aperm2',...
    'WireSec_Idensity_Aperm2','WirePriNstrands','WireSecNstrands','WirePri_per_layer',...
    'WirePri_Nlayer','WireSec_per_layer','WireSec_Nlayer','Ns1','Ns2','Ns3','Ns4',...
    'CopperPackingFactor','PackingFactor','LossCore_W','LossCopper_W','WeightCore_g',...
    'WeightPri_copper_g','WeightPri_Insu_g','WeightSec_copper_g','WeightSec_Insu_g',...
    'WeightCore_Insu_g','TotalWeight_g','TempAbsolute_C','CoreIndex'});
writetable(XfmerDesignTable,filename_xfmer,'Sheet',ResultDatasheetname);

InductorDesignTable = array2table(ResultL,'VariableNames',{'Po_W','Vin_V','Vpri_V',...
    'Vinsulation_max_V','fs_Hz','matno','CoreMatFreq_Hz','CoreCenterLegL_m','CoreCenterLegT_m',...
    'CoreAc_m2','CoreWindowH_m','CoreWindowW_m','NumOfPri','BcoreDensity_T','WirePriDia_m',...
    'WirePriFullDia_m','WirePri_Idensity_Aperm2','WirePriNstrands','WirePri_per_layer','WirePri_Nlayer',...
    'CopperPackingFactor','PackingFactor','LossCore_W','LossCopper_W','WeightCore_g','WeightPri_copper_g',...
    'WeightPri_Insu_g','WeightCore_Insu_g','TotalWeight_g','TempAbsolute_C','L','airgap_m','CoreIndex',...
    'Q','f0','A','K','RT','Ls','Cs','Cp','GT'});
writetable(InductorDesignTable,filename_inductor,'Sheet',ResultDatasheetname);

field1 = 'name';
value1_req = {'Date','Hypothesis','Notes',...
    'Q_range','f0_range','A_range','K_range',...
    'Vin_range','Vo_range','Po_range','fs_range','Winding_Pattern'};
field2 = 'data';
value2_req = {Date, Hypothesis, Notes,...
    Q_range, f0_range, A_range, K_range,...
    Vin_range, Vo_range, Po_range, fs_range, Winding_Pattern};
Requirement = struct(field1,value1_req,field2,value2_req);
Requirement_excel = squeeze(struct2cell(Requirement))';
xlswrite(filename_xfmer, Requirement_excel, Infosheetname);

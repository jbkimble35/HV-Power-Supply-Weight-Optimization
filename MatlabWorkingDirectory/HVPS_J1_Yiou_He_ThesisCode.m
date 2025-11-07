% This script is to compare different VM topology and their weight
clc, clf, clear all 
close all

% figure color
colorcodes = [0.0, 0.4470, 0.7410;...
              0.8500, 0.3250, 0.0980;...
              0.9290, 0.6940, 0.1250;...
              0,0,0];

% Save capacitance data
filename = 'Weight of different VM structure.xlsx';

% Create Structure
field1  = 'NumberOfStage';
field2  = 'TotalWeight';
field3  = 'FlyingCapacitorValue';
field4  = 'FlyingCapacitorVoltageRating';
field5  = 'NumberOfFlyingCapacitorInSeries';
field6  = 'NumberOfFlyingCapacitorInParallel';
field7  = 'OutputCapacitorValue';
field8  = 'OutputCapacitorVoltageRating';
field9  = 'NumberOfOutputCapacitorInSeries';
field10 = 'NumberOfOutputCapacitorInParallel';
field11 = 'TotalCapacitance';
field12 = 'FlyingCapEnergy';
field13 = 'OutputCapEnergy';
field14 = 'OutputImpedance';

HWCWinfo = struct(field0 ,{},field1 ,{},field2 ,{},field3 ,{},field4 ,{},field5 ,{},field6 ,{},field7,{},field8,{},field9,{},field10,{},field11,{},field12,{},field13,{});
FWCWinfo = struct(field0 ,{},field1,{},field2,{},field3,{},field4,{},field5,{},field6,{},field7,{},field8,{},field9,{},field10,{},field11,{},field12,{},field13,{});
HWDSinfo = struct(field0 ,{},field1,{},field2,{},field3,{},field4,{},field5,{},field6,{},field7,{},field8,{},field9,{},field10,{},field11,{},field12,{},field13,{});
FWDSinfo = struct(field0 ,{},field1,{},field2,{},field3,{},field4,{},field5,{},field6,{},field7,{},field8,{},field9,{},field10,{},field11,{},field12,{},field13,{});

global HWCWinfo FWCWinfo HWDSinfo FWDSinfo

% Design a voltage multiplier take input of whatever it can take to 20 kV.
Po = 375;
Vo = 20000; % half of output voltage
f  = 500000;
alpha = 0.025;
deltaVo = 50; % half of voltage ripple
Stage = 20;

for n = 1:Stage % Number of stages
    HWCWinfo(n).NumberOfStage = n;
    FWCWinfo(n).NumberOfStage = n;
    HWDSinfo(n).NumberOfStage = n;
    FWDSinfo(n).NumberOfStage = n;
    HWCWVMfactor(n) = 2*n*2;
    FWCWVMfactor(n) = n*2;
    HWDSVMfactor(n) = 2*n*2;
    FWDSVMfactor(n) = n*2;

    % Decide the capacitance based on the ripple value and voltage drop value:
    y = CapValuesWeight_updatedSSL(Po,f,Vo,n,alpha,deltaVo);
    CapHWCWWeight = y(1);
    CapHWDSWeight = y(2);
    CapFWCWWeight = y(3);
    CapFWDSWeight = y(4);
    
    if (~isempty(HWCWinfo(n).TotalCapacitance))
        HWCWCapacitance(n) = HWCWinfo(n).TotalCapacitance;
        HWCWEnergy(n) = HWCWinfo(n).FlyingCapEnergy + HWCWinfo(n).OutputCapEnergy;
        HWCWImpedance(n) = HWCWinfo(n).OutputImpedance;
    else
        HWCWCapacitance(n) = NaN;
        HWCWEnergy(n) = NaN;
        HWCWImpedance(n) = NaN;
    end
    
    if (~isempty(FWCWinfo(n).TotalCapacitance))
        FWCWCapacitance(n) = FWCWinfo(n).TotalCapacitance;
        FWCWEnergy(n) = FWCWinfo(n).FlyingCapEnergy + FWCWinfo(n).OutputCapEnergy;
        FWCWImpedance(n) = FWCWinfo(n).OutputImpedance;
    else
        FWCWCapacitance(n) = NaN;
        FWCWEnergy(n) = NaN;
        FWCWImpedance(n) = NaN;
    end
    
    if (~isempty(HWDSinfo(n).TotalCapacitance))
        HWDSCapacitance(n) = HWDSinfo(n).TotalCapacitance;
        HWDSEnergy(n) = HWDSinfo(n).FlyingCapEnergy + HWDSinfo(n).OutputCapEnergy;
        HWDSImpedance(n) = HWDSinfo(n).OutputImpedance;
    else
        HWDSCapacitance(n) = NaN;
        HWDSEnergy(n) = NaN;
        HWDSImpedance(n) = NaN;
    end
    
    if (~isempty(FWDSinfo(n).TotalCapacitance))
        FWDSCapacitance(n) = FWDSinfo(n).TotalCapacitance;
        FWDSEnergy(n) = FWDSinfo(n).FlyingCapEnergy + FWDSinfo(n).OutputCapEnergy;
        FWDSImpedance(n) = FWDSinfo(n).OutputImpedance;
    else
        FWDSCapacitance(n) = NaN;
        FWDSEnergy(n) = NaN;
        FWDSImpedance(n) = NaN;
    end
end
%% 
% Plot
% Set figure and axis positions
figure(2);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.04, 0.5, 0.85]);
ax1 = axes('Position',[0.27 0.22 0.7 0.75]);
ax1.ActivePositionProperty = 'outerposition';

figure(3);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.04, 0.5, 0.85]);
ax1 = axes('Position',[0.27 0.22 0.7 0.75]);
ax1.ActivePositionProperty = 'outerposition';

%% figure (2) plots voltage conversion ratio vs total capacitor weight of 
% the full VM at 2Po and 2Vo
FXaxis = [HWCWVMfactor; HWDSVMfactor; FWCWVMfactor; FWDSVMfactor];
FYaxis = [CapHWCWWeight; CapHWDSWeight; CapFWCWWeight; CapFWDSWeight]*2;
FXaxis(isinf(FXaxis)) = NaN;
FYaxis(isinf(FYaxis)) = NaN;
xaxislabel = 'Voltage conversion ratio ({V_{OUT}}/{V_{IN}})';
yaxislabel = 'Total weight of capacitors (kg)';
x2axislabel = 'Input voltage (kV)';
Flegend = {'Half-wave Cockcroft–Walton','Half-wave Dickson','Full-wave Cockcroft–Walton','Full-wave Dickson'};
legendx = 3;

roundXunit = 5; % label X axis every 5 unit
roundYunit = 1; % label Y axis every 10^(1) unit in log scale
figure(2);
Xplotrange = [min(min(FXaxis)), min(max(FXaxis, [], 2))];
for i = 1:1:size(FXaxis, 1)
    line(FXAxis(i,:), FYaxis(i,:), 'Color', colorcodes(i,:),'linewidth',3);hold on;
    scatter(FXAxis(i,:), FYaxis(i,:),200,colorcodes(i,:),'filled', 'MarkerFaceAlpha',1, 
        'MarkerEdgeAlpha',1); hold on;
end
hold off; grid;
set(gca, 'FontSize', 22); set(gca, 'color', 'none');

% Set x axis
xlabel(xaxislabel);
xlim(Xplotrange);
xtick = round(linspace(Xplotrange(1), Xplotrange(2), 5) / roundXunit, 0) * roundXunit; % round to the next 5
xtick(find(xtick < Xplotrange(1))) = Xplotrange(1); % cap X axis minimal
xtick(find(xtick > Xplotrange(2))) = Xplotrange(2); % cap X axis max
xticklab = strsplit(num2str(xtick));
set(gca, 'XTick', xtick, 'XTickLabel', xticklab);

% Set y axis
set(gca, 'yscale', 'log');
ylabel(yaxislabel);

YPlotIndex = intersect(find(FXaxis >= Xplotrange(1)), find(FXaxis <= Xplotrange(2)));
YPlotrange = [floor(min(min(FYaxis(YplotIndex)))/ roundYunit) * roundYunit, ...
               (floor(max(max(FYaxis(YplotIndex)))/ roundYunit) + 1) * roundYunit];
Yplotrange = [1, 4]; % uncomment for a fixed range, in this case, log for log scale only
ylim(10.^Yplotrange);
ytick = logspace(Yplotrange(1), Yplotrange(2), 4);
yticklab = strsplit(num2str(ytick));
set(gca, 'YTick', ytick, 'YTickLabel', yticklab);

% customize legend in log span
legendyspan = Yplotrange(2) - Yplotrange(1);
legendystart = legendyspan * 0.97 + Yplotrange(1);
for i = 1:1:size(Flegend,2)
    text(legendx, 10^(legendystart - 0.05*legendyspan*i), Flegend{i}, 'Color', colorcodes(i,:), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'fontsize', 22);
end
% double X axis
b = axes('Position', [0.29 0.11 .68 1e-12]);
set(b, 'Color', 'none');
set(b, 'fontsize', 22);
xlabel(b, x2axislabel)
set(b, 'xlim', [min(min(FXAxis)),min(max(FXAxis,[],2))]) % same with before
xticklab = strsplit(num2str(2 * Vo./xtick / 1000));
set(b, 'XTick', xtick, 'XTickLabel', xticklab);
%% figure (3) plots voltage conversion ratio vs total capacitor energy of 
% the full VM at 2*Po and 2*Vo
FXaxis = [HWCWVMfactor; HWDSVMfactor; FWCWVMfactor; FWDSVMfactor];
FYaxis = [HWCWEnergy; HWDSEnergy; FWCWEnergy; FWDSEnergy] / 1e12*2;
FXaxis(isinf(FXaxis)) = NaN;
FYaxis(isinf(FYaxis)) = NaN;
xaxislabel = 'Voltage conversion ratio ({V_{OUT}}/{V_{IN}})';
yaxislabel = 'Total energy in capacitors (J)';
x2axislabel = 'Input voltage (kV)';
Flegend = {'Half-wave Cockcroft–Walton','Half-wave Dickson','Full-wave Cockcroft–Walton','Full-wave Dickson'};
legendx = 3;

roundXunit = 5; % label X axis every 5 unit
roundYunit = 1; % label Y axis every 10^(1) unit in log scale
figure(3);
Xplotrange = [min(min(FXaxis)), min(max(FXaxis, [], 2))];
for i = 1:1:size(FXaxis, 1)
    line(FXAxis(i,:), FYaxis(i,:), 'Color', colorcodes(i,:),'linewidth',3);hold on;
    scatter(FXAxis(i,:), FYaxis(i,:),200,colorcodes(i,:),'filled', 'MarkerFaceAlpha',1, 
        'MarkerEdgeAlpha',1); hold on;
end
hold off; grid;
set(gca, 'FontSize', 22); set(gca, 'color', 'none');

% Set x axis
xlabel(xaxislabel);
xlim(Xplotrange);
xtick = round(linspace(Xplotrange(1), Xplotrange(2), 5) / roundXunit, 0) * roundXunit; % round to the next 5
xtick(find(xtick < Xplotrange(1))) = Xplotrange(1); % cap X axis minimal
xtick(find(xtick > Xplotrange(2))) = Xplotrange(2); % cap X axis max
xticklab = strsplit(num2str(xtick));
set(gca, 'XTick', xtick, 'XTickLabel', xticklab);

% Set y axis
set(gca, 'yscale', 'log');
ylabel(yaxislabel);
YPlotIndex = intersect(find(FXaxis >= Xplotrange(1)), find(FXaxis <= Xplotrange(2)));
YPlotrange = [floor(min(min(FYaxis(YplotIndex)))/ roundYunit) * roundYunit, ...
               (floor(max(max(FYaxis(YplotIndex)))/ roundYunit) + 1) * roundYunit];
Yplotrange = [-2, 2]; % uncomment for a fixed range, in this case, log for log scale only
ylim(10.^Yplotrange);
ytick = logspace(Yplotrange(1), Yplotrange(2), 5);
yticklab = strsplit(num2str(ytick));
set(gca, 'YTick', ytick, 'YTickLabel', yticklab);

% customize legend in log span
legendyspan = Yplotrange(2) - Yplotrange(1);
legendystart = legendyspan * 0.97 + Yplotrange(1);
for i = 1:1:size(Flegend,2)
    text(legendx, 10^(legendystart - 0.05*legendyspan*i), Flegend{i}, 'Color', colorcodes(i,:), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'fontsize', 22);
end
% double X axis
b = axes('Position', [0.27 0.11 .7 1e-12]);
set(b, 'Color', 'none');
set(b, 'fontsize', 22);
xlabel(b, x2axislabel)
set(b, 'xlim', [min(min(FXAxis)),min(max(FXAxis,[],2))]) % same with before
xticklab = strsplit(num2str(2 * Vo./xtick / 1000));
set(b, 'XTick', xtick, 'XTickLabel', xticklab);
%User inputs:
input_voltage = 400;
output_voltage = 4000;
output_power = 200;
voltage_FOS = 1.8;
current_FOS = 2;

%Topological information
num_stages = ceil(output_voltage/input_voltage);
%The number of stages in the converter
num_cap_per_stage = 3;
num_switch_per_stage = 5; %5 for Rezanejad
num_diode_per_stage = 3;
num_DCDC_per_stage = 1;
num_iso_per_stage = 3;
num_driver_per_stage = 3;

%Select lightest components possible
%Cap spreadsheet
cap_data = xlsread('Component_Masses_Capacitors');
cap_data_length = length(cap_data(:,1));
cap_masses = cap_data(:,1); %unfortunately all these are offset
cap_lengths = cap_data(:,2);
cap_widths = cap_data(:,3);
cap_voltages = cap_data(:,4);

cap_index = 1;
cap_mass = 100;
for index = 1:cap_data_length
    if cap_voltages(index) >= input_voltage*voltage_FOS
        if cap_masses(index) < cap_mass
            cap_index = index;
            cap_mass = cap_masses(index);
        end
    end
end

cap_mass = cap_masses(cap_index); %g
cap_length = cap_lengths(cap_index); %mm
cap_width = cap_widths(cap_index); %mm
cap_area = cap_length*cap_width; %mm^2

%FET spreadsheet
fet_data = xlsread('Component_Masses_FETs');
fet_data_length = length(fet_data(:,1));
fet_masses = fet_data(:,1);
fet_lengths = fet_data(:,2);
fet_widths = fet_data(:,3);
fet_voltages = fet_data(:,4);

fet_index = 1;
fet_mass = 100;
for index = 1:fet_data_length
    if fet_voltages(index) >= input_voltage*voltage_FOS
        if fet_masses(index) < fet_mass
            fet_index = index;
            fet_mass = fet_masses(index);
        end
    end
end

fet_mass = fet_masses(fet_index); %g
fet_length = fet_lengths(fet_index); %mm
fet_width = fet_widths(fet_index); %mm
fet_area = fet_length*fet_width; %mm^2

%Diode spreadsheet
diode_data = xlsread('Component_Masses_Diodes');
diode_data_length = length(diode_data(:,1));
diode_masses = diode_data(:,1);
diode_lengths = diode_data(:,2);
diode_widths = diode_data(:,3);
diode_voltages = diode_data(:,4);

diode_index = 1;
diode_mass = 100;
for index = 1:diode_data_length
    if diode_voltages(index) >= input_voltage*voltage_FOS
        if diode_masses(index) < diode_mass
            diode_index = index;
            diode_mass = diode_masses(index);
        end
    end
end

diode_mass = diode_masses(diode_index); %g
diode_length = diode_lengths(diode_index); %mm
diode_width = diode_widths(diode_index); %mm
diode_area = diode_length*diode_width; %mm^2

%DCDC spreadsheet
DCDC_data = xlsread('Component_Masses_DCDCs');
DCDC_data_length = length(DCDC_data(:,1));
DCDC_masses = DCDC_data(:,1);
DCDC_lengths = DCDC_data(:,2);
DCDC_widths = DCDC_data(:,3);
DCDC_voltages = DCDC_data(:,4);

DCDC_index = 1;
DCDC_mass = 100;
for index = 1:DCDC_data_length
    if DCDC_voltages(index) >= input_voltage*voltage_FOS
        if DCDC_masses(index) < DCDC_mass
            DCDC_index = index;
            DCDC_mass = DCDC_masses(index);
        end
    end
end

DCDC_mass = DCDC_masses(DCDC_index); %g
DCDC_length = DCDC_lengths(DCDC_index); %mm
DCDC_width = DCDC_widths(DCDC_index); %mm
DCDC_area = DCDC_length*DCDC_width; %mm^2

%iso spreadsheet
iso_data = xlsread('Component_Masses_Isos');
iso_data_length = length(iso_data(:,1));
iso_masses = iso_data(:,1);
iso_lengths = iso_data(:,2);
iso_widths = iso_data(:,3);
iso_voltages = iso_data(:,4);

iso_index = 1;
iso_mass = 100;
for index = 1:iso_data_length
    if iso_voltages(index) >= input_voltage*voltage_FOS
        if iso_masses(index) < iso_mass
            iso_index = index;
            iso_mass = iso_masses(index);
        end
    end
end

iso_mass = iso_masses(iso_index); %g
iso_length = iso_lengths(iso_index); %mm
iso_width = iso_widths(iso_index); %mm
iso_area = iso_length*iso_width; %mm^2

%driver spreadsheet
driver_data = xlsread('Component_Masses_Drivers');
driver_data_length = length(driver_data(:,1));
driver_masses = driver_data(:,1);
driver_lengths = driver_data(:,2);

driver_widths = driver_data(:,3);
driver_voltages = driver_data(:,4);

driver_index = 1;
driver_mass = 100;
for index = 1:driver_data_length
    if driver_voltages(index) >= input_voltage*voltage_FOS
        if driver_masses(index) < driver_mass
            driver_index = index;
            driver_mass = driver_masses(index);
        end
    end
end

driver_mass = driver_masses(driver_index); %g
driver_length = driver_lengths(driver_index); %mm
driver_width = driver_widths(driver_index); %mm
driver_area = driver_length*driver_width; %mm^2

%PCB material information
mass_per_area = 3.68; %kg per m^2. used www.leiten.de for mass
mass_per_mm = mass_per_area * 1e-3; %g per mm^2

%Do area calculations
raw_cap_area = num_cap_per_stage*num_stages*cap_area; %mm^2
raw_switch_area = num_switch_per_stage*num_stages*fet_area; %mm^2
raw_diode_area = num_diode_per_stage*num_stages*diode_area; %mm^2
raw_DCDC_area = num_DCDC_per_stage*num_stages*DCDC_area; %mm^2
raw_iso_area = num_iso_per_stage*num_stages*iso_area; %mm^2
raw_driver_area = num_driver_per_stage*num_stages*driver_area; %mm^2
raw_cell_area = raw_cap_area + raw_switch_area + ...
    raw_diode_area + raw_DCDC_area + raw_iso_area + raw_driver_area;
pcb_est_area = raw_cell_area*2; %gate drivers, passives, clearance

%Mass calculations
pcb_est_mass = pcb_est_area*mass_per_mm; %g
total_cap_mass = num_cap_per_stage*num_stages*cap_mass; %g
switch_mass = num_switch_per_stage*num_stages*fet_mass; %g
diode_mass = num_diode_per_stage*num_stages*diode_mass; %g

DCDC_mass = num_DCDC_per_stage*num_stages*DCDC_mass; %g
iso_mass = num_iso_per_stage*num_stages*iso_mass; %g
driver_mass = num_driver_per_stage*num_stages*driver_mass; %g
total_mass = pcb_est_mass + total_cap_mass + switch_mass + ...
    diode_mass + DCDC_mass + iso_mass + driver_mass; %g
disp(total_mass);
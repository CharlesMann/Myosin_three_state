function [time, a_attach, a_detach, Ca_conc, N_a_active] = Visualize_Actin_Rates()
% This function is used to import the actin kinetic rates from the Fortran
% threestate model. File is written as two columns. The first represents
% the a_on_rate, second is a_off_rate.

%Also going to grab the calcium

id = fopen('../Standalone_test/actin_rates.txt');
file = textscan(id, '%f%f%f');

a_attach = file{1};
a_detach = file{2};
N_a_active = file{3};
fclose(id);

id2 = fopen('../Standalone_test/calcium.txt');
file = textscan(id2, '%f%f');
Ca_conc = file{1};
time = file{2};
fclose(id2);

end
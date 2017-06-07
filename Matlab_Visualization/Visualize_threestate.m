function Visualize_threestate()
% This function will call the necessary functions to visualize text file
% data output from the Fortran standalone version of the threestate myosin
% kinetics model

% Get data
[time, a_attach, a_detach, ca_conc, N_a_active] = Visualize_Actin_Rates();
[bins, kD1_D2, kD2_D1, kD2_A, kA_D2] = Visualize_Rate_Functions();
[N_D1, N_D2, N_mbound, force] = Visualize_Pops_Props_Force();
populations = Visualize_pops();

% Plot data
% Missing Calcium
figure('Name', 'Fortran Version');

% Adding Calcium plot

subplot(8,2,[1 2]);
plot(time, ca_conc);
legend({'Ca Conc'});

% Actin on and off rates
subplot(8,2,[3 4]);
plot(time, a_attach, 'b');
hold on;
plot(time, a_detach, '--r');
legend({'Actin on', 'Actin off'});


% Activated sites and Number bound
% Need to add activated sites
subplot(8,2,[5 6]);
plot(time, N_mbound, 'r');
hold on;
plot(time, N_a_active, 'b');
ylim([0 1]);
legend({'N_mbound', 'Actin Active'});



%Bin Popoulation distribution
subplot(8,2,[7 8]);
for i = 1:numel(time)
    if mod(i,100) ==1
        cla
        plot(bins, populations(i,:));
        legend({'Bin Distributions'});
        ylim([0 1]);
        drawnow;
    end
end

% D2 -> A and A-> D2 rates
subplot(8,2,[9 10]);
for i = 1:numel(time)
    if mod(i,100) == 1
        cla
        plot(bins, kD2_A(i,:), 'b');
        hold on;
        plot(bins, kA_D2(i,:), 'r');
        legend({'D2->A','A->D2'});
        ylim([0 20]);
    end
end


% D1 and D2 rates
% Double check that the D1/D2 rates are correct
subplot(8,2,[11 12]);
plot(time, kD1_D2, 'b');
hold on;
plot(time, kD2_D1, 'r');
legend({'D1->D2', 'D2->D1'});



% D1 and D2 Populations
subplot(8,2,[13 14]);
plot(time, N_D1, 'b');
hold on;
plot(time, N_D2, 'r');
legend({'D1', 'D2'});


% HS Force
subplot(8,2,[15 16]);
plot(time, force);
legend({'HS Force'});






end
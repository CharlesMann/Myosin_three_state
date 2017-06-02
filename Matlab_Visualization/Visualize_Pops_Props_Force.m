function [N_D1, N_D2, N_mbound, force] = Visualize_Pops_Props_Force()
% View populations in D1, D2, the number bound, and the force generated at
%each time step. Need population distribution by bins

id = fopen('../Standalone_test/populations_force.txt');
file = textscan(id, '%f%f%f');
N_D1 = file{1};
% Column 1 is D1, 2 is D2, 3 is N_m_bound, 4 is force
%N_D1 = file{1};
N_D2 = file{2};
force = file{3};

id2 = fopen('../Standalone_test/mbound.txt');
file2 = textscan(id2,'%f');
N_mbound = file2{1};

fclose(id);
% Will add bin populations

end
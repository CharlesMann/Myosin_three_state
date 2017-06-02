function [populations] = Visualize_pops()

id = fopen('../Standalone_test/pop_distribution.txt');
pops = cell2mat(textscan(id, '%f'));
for i = 1:size(pops)/41
    populations(i,:) = pops((i-1)*41+1:i*41,1);
end


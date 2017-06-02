function [bins, kD1_D2, kD2_D1, kD2_A, kA_D2] =  Visualize_Rate_Functions()
% Input rate functions k1-k4, plot with bin displacements

bins = -10:.5:10;
id1 = fopen('../Visual Studio Debug Attempt/kD1_D2.txt');
kD1_D2 = cell2mat(textscan(id1, '%f'));

id2 = fopen('../Visual Studio Debug Attempt/kD2_D1.txt');
kD2_D1 = cell2mat(textscan(id2, '%f'));

id3 = fopen('../Visual Studio Debug Attempt/kD2_A.txt');
kD2_A_raw = cell2mat(textscan(id3, '%f'));

% This now needs to be separated into chunks of 41
kD2_A = zeros(numel(kD2_A_raw)/41,41);
for i = 1:size(kD2_A,1)
    kD2_A(i,:) = kD2_A_raw((i-1)*41+1:i*41);
end


id4 = fopen('../Visual Studio Debug Attempt/kA_D2.txt');
kA_D2_raw = cell2mat(textscan(id4, '%f'));
kA_D2 = zeros(numel(kA_D2_raw)/41,41);
for i = 1:size(kA_D2,1)
    kA_D2(i,:) = kA_D2_raw((i-1)*41+1:i*41);
end

fclose(id1);
fclose(id2);
fclose(id3);
fclose(id4);

end
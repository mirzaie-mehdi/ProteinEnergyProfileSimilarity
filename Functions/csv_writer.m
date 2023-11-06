clear all
clc
load('H:\Dunbrack\new_energy_dis5\energy_dell_Dis5.mat')
fid =  fopen('H:\Dunbrack\new_energy_dis5\energy_dell_Dis5.csv','wt');
for i =1: 167
    disp(num2str(i))
    %[atnamei,amnamei]=num167toatomtype(i);
    for j = i:167
%         [atnamej,amnamej]=num167toatomtype(j);
%         fprintf(fid,'%s \t', amnamei);
%         fprintf(fid,'%s \t', atnamei);
%         fprintf(fid,'%s \t', amnamej);
%         fprintf(fid,'%s \t', atnamej);
        for k = 1:9
            fprintf(fid,'%20.15f \t', energy_dell_Dis5{i,j}(1,k));
        end
        fprintf(fid,'\n');
    end
end
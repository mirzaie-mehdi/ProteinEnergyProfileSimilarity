clear;clc
load('zzi.mat');
Pij = zeros(20,20);
for i=1:20
    disp(num2str(i))
    v = variable(i);
    grd = gradient(zzi{1,i},v);
    slv = solve(grd);
    slv_cel= struct2cell(slv);
    slv_cel2=slv_cel([5 10 8 18 19 20 11 2 7 17 15 1 12 14 13 16 3 6 9 4]);
    for j=1:20
        Pij(i,j) = vpa(slv_cel2{j},10);
    end
end
writematrix(Pij,'Pij.csv')

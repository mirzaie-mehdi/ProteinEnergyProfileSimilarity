clear;clc
path='../Data/Dunbrack/';
load([path 'zzi_torsion2.mat']);
Pi = zeros(1,20);
for i=1:20
    disp(num2str(i))
    v = variable2(i);
    grd = gradient(zzi{1,i},v);
    slv = solve(grd);
    Pi(1,i) = vpa(slv,10);
end
writematrix(Pi,[path '/Pi_torsion2.csv'])
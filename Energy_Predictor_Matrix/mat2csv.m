load('../Train_Energy/energy_dell_dunbrack.mat')
rows = 167;
cols = 167;
numColumns = 29;
numRows = rows * (rows + 1) / 2;
resultMatrix = zeros(numRows, numColumns);
rowCounter = 1;
for i = 1:rows
    disp(num2str(i))
    for j = i:cols
        resultMatrix(rowCounter, :) = energy_dell_dunbrack{i, j};
        rowCounter = rowCounter + 1;
    end
end
writematrix(resultMatrix,'energy_dell_dunbrack_rm_Olaps.csv')
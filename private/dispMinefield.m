% Display entire minefield in figure window
function dispMinefield()
    global minefieldDim;

    for i = 1:minefieldDim(1)
        for j = 1:minefieldDim(2)
            dispCell(i,j);
        end
    end
end
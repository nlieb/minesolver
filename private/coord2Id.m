%Converts a coordinate on the minefield to a position ID
function id = coord2Id(row, col)
    global minefieldDim
    id = (row-1) * minefieldDim(2) + col;
end
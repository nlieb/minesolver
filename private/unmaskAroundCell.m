%Unmasks around a cell
function unmaskAroundCell(row, col)
    global minefield minefieldDim
    
    for m = (row-1):(row+1)
        for n = (col-1):(col+1)
            if(m > 0 && n > 0 && m <= minefieldDim(1) && n <= minefieldDim(2))
                minefield(m,n,3) = 0;
            end
        end
    end
end
%Returns the number of mines currently solved in the minefield
function mines = minesSolved()
    global minefield minefieldDim
    
    mines = 0;
    
    for i = 1:minefieldDim(1)
        for j = 1:minefieldDim(2)
            if minefield(i,j,2) == 99
                mines = mines+1;
            end
        end
    end
end
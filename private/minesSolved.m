%Returns the number of mines currently solved in the minefield
function nMines = minesSolved()
    global minefield
    
    %Vectorized approach to finding marked mines
    mines = find(minefield(:,:,2) == 99);
    nMines = length(mines);
end
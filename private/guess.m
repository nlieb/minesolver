function success = guess()
    %Guesses the cell with the lowest probability of being a mine
    global minefield minefieldDim mineNum
    global probVector
    
    fprintf('Initiating probability calculations...\n');
    
    %Find the number of bombs solved
    bombsSolved = minesSolved();
    
    %Create a probability matrix
    probVector = zeros(1, minefieldDim(1)*minefieldDim(2));
    
    %Fill the probability matrix for unknown cells touching hints
    for pos = 1:length(probVector)
        coord = id2Coord(pos);
        if minefield(coord(1),coord(2),2) ~= -1 && minefield(coord(1),coord(2),2) ~= 99
            probability = calcProbability(coord(1),coord(2));
            if probability ~= 0
                addProbability(coord(1),coord(2),probability)
            end
        end
    end
    
    voidVector = [];
    %Fill the probability matrix for cells not touching hints
    for pos = 1:length(probVector)
        coord = id2Coord(pos);
        if minefield(coord(1),coord(2),2) == -1 && probVector(pos) == 0
            voidVector = [voidVector, pos];
        end
    end
    
    %This is the probability of guessing a mine in the "void", i.e not touching a clue
    voidProbability = length(voidVector)/(mineNum-bombsSolved);
    
    %Apply this to the probability vector
    for j = 1:length(voidVector)
        probVector(voidVector(j)) = voidProbability;
    end
    
    %We must first replace all 0s by Inf to find the non-zero minimum element
    probVector(probVector == 0) = Inf;
    
    %Pick a cell to guess
    low = find(probVector == min(probVector));
%     pick = ceil(rand*length(low));
%     cellCoord = id2Coord(low(pick)); 
    cellCoord = id2Coord(low(1));
    successChance = (1 - min(probVector))*100;
    
    fprintf('Guessing at (%2d,%2d)\nSuccess probability: %4.1f%%\n',cellCoord(1),cellCoord(2),successChance);
    
    %Unmask the guess
    minefield(cellCoord(1),cellCoord(2),3) = 0;
    updateMaskedMinefield();
    
    %Check if the guess is a mine (return if true)
    if minefield(cellCoord(1),cellCoord(2),2) == 9
        fprintf('Oh no! I guessed wrong!\n');
        success = false;
    else
        %Guess is successful, try to continue solving
        success = true;
    end
end

function addProbability(row, col, probability)
    %Appends a probability around a cell
    global minefield minefieldDim
    global probVector
    
    for m=(row-1):(row+1)
        for n=(col-1):(col+1)
            %Check bounds
            if m <= minefieldDim(1) && n <= minefieldDim(2) && m > 0 && n > 0 && ~(m == row && n == col)
                if minefield(m,n,2) == -1
                    pos = coord2Id(m,n);
                    %Replace only with larger probability
                    if probability > probVector(pos)
                        probVector(pos) = probability;
                    end
                end
            end
        end
    end
end

function probability = calcProbability(row, col)
    %Calculates the probability of there being a mine around a cell
    global minefield minefieldDim
    
    mineCounter = 0;
    unknownCounter = 0;
    
    %Cycle around to count mines and unknowns
    for m=(row-1):(row+1)
        for n=(col-1):(col+1)
            if m <= minefieldDim(1) && n <= minefieldDim(2) && m > 0 && n > 0
                if minefield(m,n,2) == 99
                    mineCounter = mineCounter+1;
                elseif minefield(m,n,2) == -1
                    unknownCounter = unknownCounter+1;
                end
            end
        end
    end
    
    %Default probability is zero
    probability = 0;
    
    %Must have unknown around to calculate probability
    if unknownCounter == 0
        return;
    else
        probability = (double(minefield(row,col,2)) - mineCounter)/unknownCounter;
    end
end



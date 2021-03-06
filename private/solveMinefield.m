% Contains main solving loop.
function solveMinefield()
    global minefieldDim mineNum
    global equationMatrix equationMatrixPos
    global bombs equations
    global grain solvedArray
    
    %Initialize variable used to check if the matrix is solved (i.e no
    %changes have been made)
    equations = 1;
    lastPassBombs = -1;
    bombs = 0;
    
    while (equations > 0 || lastPassBombs ~= bombs)
        %Reset solved checks
        lastPassBombs = bombs;
        bombs = 0;
        equations = 0;
        equationMatrixPos = 1;
        
        %Run a first pass of the minefield with the basic expansion algorithm
        if(clearPass() ~= 0)
            equations = equations+1;
        end
        
        %Build the equationMatrix
        for m = 1:minefieldDim(1)
            n = 0;
            while(n < minefieldDim(2))
                n = n+1;
                if(mod(n-1, grain) == 0 && solvedArray(ceil(m/grain),ceil(n/grain)))%dnc optimization
                    n = n+grain;%skip to next block
                    continue;
                end
                getEquationBuilder(m,n);
            end
        end
        
        %Solve and parse the equationMatrix and SolvedEqMatrix
        solveEquations(lastPassBombs);
        
        %Implement solved changes to unmask minefield
        updateMaskedMinefield();
        
        %Zero out the equation matrix
        equationMatrix(:,:) = 0;
        
        %Run mine counter if nothing was solved yet
        if (equations == 0 && lastPassBombs == bombs)
            if(bombs ~= mineNum)
                fprintf('running mine counting method\n');
                
                mineCountingMethod(lastPassBombs);
                updateMaskedMinefield();
                
                %Zero out the equation matrix
                equationMatrix(:,:) = 0;
            end
        end
        
        %Display masked minefield
        %disp(minefield(:,:,2));
        drawnow;
        dncRecheckMinefield();
        fprintf('bombs found: %d\n', bombs);
    end
    
    %fprintf('\nDeterministic solving complete!\n\n');
end

%Check to see if any squares can be solved outright
function counter = clearPass()
    global minefield minefieldDim 
    global grain solvedArray

    isUpdated = true;
    counter = 0;
    
    drawnow;
    while(isUpdated)
        if(mod(counter,2) ~= 0)
            drawnow;
        end
        dncRecheckMinefield(); %refresh dnc array
        
        temp = minefield(:,:,2);
        
        for m = 1:minefieldDim(1)
            n = 1;
            while n < minefieldDim(2)
                if(mod(n-1, grain) == 0 && solvedArray(ceil(m/grain),ceil(n/grain)))%dnc optimization
                    n = n+grain;%skip to next block
                    continue;
                end
                if(minefield(m,n,2) ~= -1 && minefield(m,n,2) ~= 99)
                    simpleSolve(m,n);
                end
                n = n+1;
            end
        end

        updateMaskedMinefield();
        
        %Check if any changes were made
        if(temp == minefield(:,:,2))
            isUpdated = false;
        else
            counter = counter + 1;
        end
    end
    drawnow;
end

%Deals with the cases where a cell's unknowns are ALL mines or ALL safe
function simpleSolve(row, col)
    global minefield minefieldDim
    
    %Shortcut 0s to clear
    if(minefield(row,col,2) == 0)
        %status = 'All safe';
        unmaskAroundCell(row,col);
        return;
    end
  
    mines = 0;
    unknowns = 0;
    
    %Count the number of unknowns and mines around the cell
    for m = (row-1):(row+1)
        for n = (col-1):(col+1)
            if(m > 0 && n > 0 && m <= minefieldDim(1) && n <= minefieldDim(2))
                if(minefield(m,n,2) == -1)
                    unknowns = unknowns + 1;
                elseif(minefield(m,n,2) == 99)
                    mines = mines + 1;
                end
            end
        end
    end
    
    %Check if the number of mines is equal to the hint
    if(mines > 0 && minefield(row,col,2) == mines)
        %status = 'All safe'
        unmaskAroundCell(row,col);
        return;
    end
    
    %Check if all the unknowns are mines
    if((minefield(row,col,2) - unknowns - mines) == 0)
        %status = 'All mines'
        buildSolvedEquationsAroundCell(row,col);
        return;
    end
    
    %status = 'Undefined';
end

%Builds then solves the solver matrix
function solveEquations(lastBombs)
    global mineNum
    global solvedEqMatrix equationMatrix
    global solverMatrix sizeS
    global equations bombs
    global maxPermutations

    %Concatonate the solved and equation matrices vertically to set up
    %the row reduction
    solverMatrix = vertcat(solvedEqMatrix, equationMatrix);
    sizeS = size(solverMatrix);
    
    %dispEquations();

    %Row reduce the matrix (i.e solve it)
    solverMatrix = frref(solverMatrix);
    %Round the results to prevent rounding errors
    solverMatrix = round(solverMatrix);
    
    %dispEquations();
    
    %Parse equation Matrix for solved rows
    parseEquations();
    
    iteration = 0;
    
    while equations == 0 && bombs == lastBombs && iteration < maxPermutations && bombs ~= mineNum
        equations = 0;
        bombs = 0;
        iteration = iteration + 1;
        
        if iteration == 1
            permVector = 1:sizeS(2);
            %Create and store the permuted matrix
            permMatrix = zeros(sizeS(1),sizeS(2));
        end
        
        %Permute the columns of the solvermatrix
        if(mod(iteration,2) ~= 0)
           permVector = [permVector(sizeS(2)-1:-1:1),sizeS(2)];%flip the permute vector
        else
           permVector = [randperm(sizeS(2)-1), sizeS(2)];%randomly permute the columns
        end

        for j = 1:sizeS(2)
            permMatrix(:,j) = solverMatrix(:,permVector(j));
        end
        
        %rref the matrix
        permMatrix = frref(permMatrix);
        %Round the results to prevent rounding errors
        permMatrix = round(permMatrix);
        
        %Convert back to regular id order in solverMatrix
        for j = 1:sizeS(2)
            solverMatrix(:,permVector(j)) = permMatrix(:,j);
        end
        
        %dispEquations();
        
        %Parse the solverMatrix
        parseEquations();
    end
    
    if(iteration>0)
        fprintf('permuted %d times\n', iteration);
    end
end

%Marks a mine in the solved equation matrix (i.e with one variable = 1)
function buildSolvedEquation(row, col)
    global solvedEqMatrix solvedEqMatrixDim
    global minefield
    
    %If yes, add solved equation to solvedEqMatrix
    id = coord2Id(row,col);
    solvedEqMatrix(id, id) = 1;
    solvedEqMatrix(id, solvedEqMatrixDim(2)) = 1;
    
    minefield(row, col, 2) = 99; %Mark cell as a mine
end

%A subfunction that builds solved equations at all unknowns around a mine
function buildSolvedEquationsAroundCell(row, col)
    global minefield minefieldDim;
    
    for m = (row-1):(row+1)
        for n = (col-1):(col+1)
            if(m <= minefieldDim(1) && n <= minefieldDim(2) && m > 0 && n > 0 && ~(m == row && n == col))
                if(minefield(m,n,2) == -1)
                    buildSolvedEquation(m,n);
                end
            end
        end
    end
end

%Checks if an equation should be built at a certain cell (i.e is touching
%at least one unknown)
function getEquationBuilder(row, col)
    global minefield minefieldDim
    
    if(minefield(row,col,2) ~= -1 && minefield(row,col,2) ~= 99)
        %Cycle around the queried cell
        for m = (row-1):(row+1)
            for n = (col-1):(col+1)
                %Check bounds
                if(m <= minefieldDim(1) && n <= minefieldDim(2) && m > 0 && n > 0 && ~(m == row && n == col))
                    %Does cell have an unknown or mine around it?
                    if(minefield(m,n,2) == -1)
                        %If yes, add equation to equationMatrix and exit function
                        buildEquation(row,col);
                        return;
                    end
                end
            end 
        end
    end
end

%Inserts information into the equation matrix
function buildEquation(row, col)
    global minefield minefieldDim
    global equationMatrix equationMatrixPos equationMatrixDim
    
    %Cycle around cell
    for m=(row-1):(row+1)
        for n=(col-1):(col+1)
            %Check bounds and if it is an unknown or mine
            if(m <= minefieldDim(1) && n <= minefieldDim(2) && m > 0 && n > 0 && (minefield(m,n,2) == -1 || minefield(m,n,2) == 99) && ~(m == row && n == col))
                %Set linear equation coefficient to one 
                equationMatrix(equationMatrixPos, coord2Id(m,n)) = 1;
            end
        end
    end
    
    equationMatrix(equationMatrixPos, equationMatrixDim(2)) = minefield(row,col,2);
    equationMatrixPos = equationMatrixPos+1;
end

%Solves minefield with added equation where all unknowns = mines left
function mineCountingMethod(lastPassBombs)
    global minefield minefieldDim mineNum
    global equationMatrix equationMatrixPos equationMatrixDim
    global bombs
    global grain solvedArray
    
    bombsLeft = mineNum - minesSolved();
    bombs = 0;
    equationMatrixPos = 1;
                
    %Build the equationMatrix
    for m = 1:minefieldDim(1)
        n = 0;
        while(n < minefieldDim(2))
            n = n+1;
            if(mod(n-1, grain) == 0 && solvedArray(ceil(m/grain),ceil(n/grain)))%dnc optimization
                n = n+grain;%skip to next block
                continue;
            end
            getEquationBuilder(m,n);
        end
    end
    
    for m = 1:minefieldDim(1)
        for n = 1:minefieldDim(2)
            %Check if it is an unknown
            if minefield(m,n,2) == -1
                %Set linear equation coefficient to one 
                equationMatrix(equationMatrixPos, coord2Id(m,n)) = 1;
            end
        end
    end
    
    %The final value is equal to the number of mines not yet found
    equationMatrix(equationMatrixPos, equationMatrixDim(2)) = bombsLeft;
    equationMatrixPos = equationMatrixPos+1;
    
    %dispEquations();
    solveEquations(lastPassBombs);
    %dispEquations();
end

%Parse the solverMatrix for solved rows
%calls the two methods, allClearMethod, and allMatchMethod
function parseEquations()
    global solverMatrix sizeS 
    global equations
    
    for i = 1:sizeS(1)
        %case when all variables don't have mines, we are
        %looking for equations of the form 1...1...|0 or -1...-1...|0
        if solverMatrix(i, sizeS(2)) == 0
            equations = equations + allClearMethod(i);
            continue; %skips the other cases below
        end
        
        %case when terminator isn't 0, we are
        %looking for equations of the form 1...1...|2 or -1...-2...|-3 etc.
        if abs(solverMatrix(i, sizeS(2))) ~= 0
            allMatchMethod(i);
        end
    end
end

%special case when all variables don't have mines, we are looking for equations of the form 1...1...|0 or -1...-1...|0
function equations = allClearMethod(i)
    global minefield
    global solverMatrix sizeS
    
    j = 1;
    equations = 0;
    
    %Loop through until you find a non-zero term
    while(j < sizeS(2) && solverMatrix(i, j) == 0)
        j = j + 1;
    end
    if j == sizeS(2) %If all terms in that row were zero, exit
        return;
    end
    check = solverMatrix(i, j);
    
    %Check if all the remaining variables are either 0 or
    %the same coefficient as the check
    n = length(find(solverMatrix(i,:) ~= check & solverMatrix(i,:) ~= 0));
    if n %n must be 0 to satisfy the form.
        return;
    end
    
    solved = find(solverMatrix(i,:) == check);
    for k = 1:length(solved)%Unmask the cell
        coord = id2Coord(solved(k));
        minefield(coord(1),coord(2),3) = 0; 
    end
    equations = equations + length(solved);
    
    return 
end

%Checks and solves equations of the form 1...1...|2 or 1...1...1|3 etc.
function allMatchMethod(i)
    global minefield
    global solverMatrix sizeS
    global equations bombs
    
    total = solverMatrix(i,sizeS(2));    
    
    if(total > 0)
        matches = find(solverMatrix(i,1:sizeS(2)-1) > 0);
        plus = true;
    else
        matches = find(solverMatrix(i,1:sizeS(2)-1) < 0);
        plus = false;
    end

    %Add up the sum of the matches
    counter = sum(solverMatrix(i,matches));
    
    %Unmask/Mark cells
    if(counter == total)
        %Look for non-zero non-matches in the solverMatrix
        if plus
            nMatches = find(solverMatrix(i,1:sizeS(2)-1) < 0);
        else
            nMatches = find(solverMatrix(i,1:sizeS(2)-1) > 0);
        end      
        
        for index = 1:length(matches)
            coord = id2Coord(matches(index));
            buildSolvedEquation(coord(1),coord(2)); %Build an equation for the mine
            bombs = bombs+1;
        end
        for index = 1:length(nMatches)
            coord = id2Coord(nMatches(index));
            equations = equations + 1;
            minefield(coord(1),coord(2),3) = 0; %Unmask the cell
        end
    end
end



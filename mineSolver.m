%Randomly solves a minesweeper field programatically.
function mineSolver
    clear; clc; close all;
    addpath('lib/', 'img/');
    
    global minefield minefieldDim mineNum
    global equationMatrix equationMatrixDim equationMatrixPos
    
    %Set minefield dimensions
    minefieldDim(1) = 60;
    minefieldDim(2) = 110;
    mineNum = 1000;
    
    initializeFigureWindow(0.2);
    dncInit();
    
    %generate the minefield into global variable "minefield"
    generateMinefield(minefieldDim(1), minefieldDim(2), mineNum);
    
    
    %save('field.mat');
    %load('field.mat');
    
    
    % Get total number of cells
    cellNum = minefieldDim(1)*minefieldDim(2);
    
    %Create global matrix to store calculations
    %Size dynamically grows
    equationMatrix = zeros(1,cellNum+1);
    equationMatrixDim = size(equationMatrix);
    equationMatrixPos = 1;
    
    solveMinefield();
end

%initialize divide and conquor data structures
function dncInit()
    global minefieldDim solvedArray grain;
    
    grain = 3;
    saM = ceil(minefieldDim(1)/grain);
    saN = ceil(minefieldDim(2)/grain);
    
    %solvedArray is the solved array which represents each 8*8 block as a flag
    %denoting whether the block is solved or not.
    solvedArray = zeros(saM, saN);
end

%Initialize the figure window with given size, offset and title
function initializeFigureWindow(multiplier)
    %Creates a figure window with defined properties
    global minefieldDim sprites
    
    sprites = {};
    sprites{1} = imread('OneCell.png');
    sprites{2} = imread('TwoCell.png');
    sprites{3} = imread('ThreeCell.png');
    sprites{4} = imread('FourCell.png');
    sprites{5} = imread('FiveCell.png');
    sprites{6} = imread('SixCell.png');
    sprites{7} = imread('SevenCell.png');
    sprites{8} = imread('EightCell.png');
    sprites{9} = imread('MineCell.png');
    sprites{10} = imread('UnknownCell.png');
    sprites{11} = imread('ZeroCell.png');
    sprites{12} = imread('FlagCell.png');
    
    %Get the screensize and figure size
    screenSize = get(0, 'ScreenSize');
    %44x44 for full size
    figureSize = [minefieldDim(2) * multiplier*44, minefieldDim(1) * multiplier*44];
    
    %Get the centered offset for the figure window
    xpos = (screenSize(3)-figureSize(1))/2;
    ypos = (screenSize(4)-figureSize(2))/2;
    
    figure('Position', [xpos, ypos, figureSize(1), figureSize(2)], ...
           'Name', 'MineSolver', 'NumberTitle', 'off', ...
           'Resize', 'off', ...
           'Toolbar', 'none', ...
           'Menubar', 'none');
       
    hold on;
    axis([0, minefieldDim(2)*35, 0, minefieldDim(1)*35]);
    axis('off');
end

% Contains main solving loop.
function solveMinefield()
    global minefield minefieldDim mineNum
    global equationMatrix equationMatrixDim equationMatrixPos
    global solvedEqMatrix solvedEqMatrixDim
    global bombs equations
    global solvedArray

    %Create the solved equation matrix
    solvedEqMatrix = zeros(equationMatrixDim(2),equationMatrixDim(2));
    solvedEqMatrixDim = size(solvedEqMatrix);
    
    %Initialize variable used to check if the matrix is solved (i.e no
    %changes have been made)
    equations = 1;
    lastPassBombs = -1;
    bombs = 0;
    
    %Run a first pass of the minefield with the basic expansion algorithm
    dispMinefield();
    clearPass();
    dispMinefield();
    
    while (equations > 0 || lastPassBombs ~= bombs)
        %Reset solved checks
        lastPassBombs = bombs;
        bombs = 0;
        equations = 0;
        equationMatrixPos = 1;
        
        if(clearPass() ~= 0)
            equations = equations+1;
        end
        
        %Build the equationMatrix
        for m = 1:minefieldDim(1)
            for n = 1:minefieldDim(2)
                getEquationBuilder(m,n);
            end
        end
        
        %Solve and parse the equationMatrix and SolvedEqMatrix
        solveEquations();
        
        %Implement solved changes to unmask minefield
        updateMaskedMinefield();
        
        %Zero out the equation matrix
        equationMatrix(:,:) = 0;
        
        %Run mine counter if nothing was solved yet
        if (equations == 0 && lastPassBombs == bombs)
            if(bombs ~= mineNum)
                fprintf('Running mine counting method:\n');
                bombsLeft = bombs;
                equationMatrixPos = 1;
                bombs = 0;
                
                mineCountingMethod(bombsLeft);
                
                %Build the equationMatrix
                for m = 1:minefieldDim(1)
                    for n = 1:minefieldDim(2)
                        getEquationBuilder(m,n);
                    end
                end
        
                solveEquations();
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
    
    fprintf('\n Deterministic solving complete!\n');
end

%Check to see if any squares can be solved outright
function counter = clearPass()
    global minefield minefieldDim solvedArray grain;

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
            while n <= minefieldDim(2)
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
function solveEquations()
    global mineNum
    global solvedEqMatrix
    global equationMatrix
    global bombs equations
    global lastBombs

    %Concatonate the solved and equation matrices vertically to set up
    %the row reduction
    solverMatrix = vertcat(solvedEqMatrix, equationMatrix);
    sizeS = size(solverMatrix);
    
    %displayEquations(solverMatrix, size(solverMatrix));

    %Row reduce the matrix (i.e solve it)
    solverMatrix = frref(solverMatrix);

    %displayEquations(solverMatrix, size(solverMatrix));
    
    %Parse equation Matrix for solved rows
    parseEquations(solverMatrix, sizeS);
    
    iteration = 0;
    
    while equations == 0 && bombs == lastBombs && iteration < 10 && bombs ~= mineNum
        equations = 0;
        bombs = 0;
        iteration = iteration + 1;
        
        %Permute the solvermatrix
        permVector = [randperm(sizeS(2)-1), sizeS(2)];
        
        %Create and store the permuted matrix
        permMatrix = zeros(sizeS(1),sizeS(2));

        for j = 1:sizeS(2)
            permMatrix(:,j) = solverMatrix(:,permVector(j));
        end
        
        %rref the matrix
        permMatrix = frref(permMatrix);
        
        %Convert back to regular id order in solverMatrix
        for j = 1:sizeS(2)
            solverMatrix(:,permVector(j)) = permMatrix(:,j);
        end
        
        %Parse the solverMatrix
        parseEquations(solverMatrix, sizeS);
    end
    
    lastBombs = bombs;
end

%Parse the solverMatrix for solved rows
function parseEquations(solverMatrix, sizeS)
    %Parse equation Matrix for solved rows
    global minefield minefieldDim
    global equations
    
    solvedSquares = zeros(minefieldDim(1)*minefieldDim(2));
    
    for i = 1:sizeS(1)

        %special case when all variables don't have mines, we are
        %looking for equations of the form 1...1...|0 or -1...-1...|0
        if(solverMatrix(i, sizeS(2)) == 0)
            j = 1;
            exit = false;

            %Loop through until you find a non-zero term
            while(j < sizeS(2) && solverMatrix(i, j) == 0)
                j = j + 1;
            end

            if j == sizeS(2) %If all terms in that row were zero, exit
                exit = true;
            else
                %Check if all the remaining variables are either 0 or
                %the same coefficient
                check = solverMatrix(i, j);

                %Loop through the rest of the array
                while j < sizeS(2)
                    if solverMatrix(i, j) ~= check && solverMatrix(i, j) ~= 0
                        exit = true;
                        break;
                    end
                    j = j + 1;
                end
            end 

            if(~exit)
                for j = 1:sizeS(2)-1
                    if solverMatrix(i, j) ~= 0
                        solvedSquares(j) = 1;

                        %coord = id2Coord(j);
                        %fprintf('(%d, %d): is empty (constant = 0)\n', coord(1), coord(2));
                    end
                end

                for k = 1:minefieldDim(1)*minefieldDim(2)% unmask special case when all variables are don't have mines
                    if(solvedSquares(k) == 1)
                        coord = id2Coord(k);
                        equations = equations + 1;
                        minefield(coord(1),coord(2),3) = 0; %Unmask the cell
                    end
                end
            end

            continue; %skips the other cases below
        end
        
        %special case when terminator isn't 0, we are
        %looking for equations of the form 1...1...|2 or -1...-1...-1|-3 etc.
        
        if(abs(solverMatrix(i, sizeS(2))) ~= 0)
            allMatchMethod(i, solverMatrix, sizeS);
        end
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
function mineCountingMethod(bombsLeft)
    global minefield minefieldDim mineNum
    global equationMatrix equationMatrixPos equationMatrixDim
    
    for m = 1:minefieldDim(1)
        for n = 1:minefieldDim(2)
            %Check bounds and if it is an unknown or mine
            if minefield(m,n,2) == -1
                %Set linear equation coefficient to one 
                equationMatrix(equationMatrixPos, coord2Id(m,n)) = 1;
            end
        end
    end
    
    %The final value is equal to the number of mines not yet found
    equationMatrix(equationMatrixPos, equationMatrixDim(2)) = mineNum - bombsLeft;
    equationMatrixPos = equationMatrixPos+1;
end

%Checks and solves equations of the form 1...1...|2 or 1...1...1|3 etc.
function allMatchMethod(i, solverMatrix, sizeS)
    global minefield
    global equations bombs
    
    total = solverMatrix(i,sizeS(2));
    counter = 0;
    
    
    if(total > 0)
        match = 1;
    else
        match = -1;
    end
    
    %Loop through the solverMatrix
    for j = 1:sizeS(2)-1
        if (solverMatrix(i, j) == match)
            counter = counter+1;
        end
    end
    
    %Unmask/Mark cells
    if(counter == total)
        for j = 1:sizeS(2)-1
            if solverMatrix(i, j) == match
                coord = id2Coord(j);
                buildSolvedEquation(coord(1),coord(2));
                bombs = bombs+1;
            elseif solverMatrix(i, j) ~= 0
                coord = id2Coord(j);
                equations = equations + 1;
                minefield(coord(1),coord(2),3) = 0; %Unmask the cell
            end
        end
    end
end

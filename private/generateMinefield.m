%Generates the minefield

function generateMinefield(m, n, mines)
    global minefield minefieldDim mineNum
    
    %Create the complete minefield on page 1
    minefield = zeros(m,n,'int8');
    minefieldDim = size(minefield);
    mineNum = mines;
    
    %Create a second page to contain the masked minefield
    minefield(:,:,2) = ones(m,n,'int8')*-1;
    
    %Apply masking on a third page
    %This will keep track of which values are shown to the AI 
    minefield(:,:,3) = ones(m,n,'int8');
    
    %Set point in the center of the minefield that will be clear and unmasked 
    startPoint = [ceil(m/2), ceil(n/2)];
    %Unmask starting point and cells around it
    minefield(startPoint(1),startPoint(2),3) = 0;
    unmaskAroundCell(startPoint(1), startPoint(2));
    
    isGen = 0;
    %Generate mine distribution
    while(isGen < mines)
        %Pick random location to place mine
        row = ceil(rand*m);
        col = ceil(rand*n);
        
        %Check if a mine isn't already present at that location
        if (minefield(row, col, 1) ~= 9)
            %Don't place mines at or around the start point (default at
            %center)
            if(row > startPoint(1)+1 || row < startPoint(1)-1 || ...
               col > startPoint(2)+1 || col < startPoint(2)-1)
                %Place mine
                minefield(row, col, 1) = 9;
                %increment the counter
                isGen = isGen+1;
            end
        end
    end
    
    %Generate the clues
    addClues();
    updateMaskedMinefield();
end

%Adds the clues to the minefield
function addClues()
   global minefield
   
   [row, col] = find(minefield(:, :, 1) == 9);
   
   for i=1:length(row);
       %Increment clues around each mine +1
       addMine(row(i), col(i));
   end
end

%Increments a hint in the minefield
function addMine(row, col)
    global minefield minefieldDim
    
    %Cycle around mine
    for m = (row-1):(row+1)
        for n = (col-1):(col+1)
            %Check bounds
            if(m > 0 && n > 0 && m <= minefieldDim(1) && n <= minefieldDim(2)) 
                if(minefield(m,n,1) ~= 9)
                    %Add mine
                    minefield(m,n,1) = minefield(m,n,1) + 1;
                end
            end
        end
    end
end
%Updates the masked minefield with the latest masking
function updateMaskedMinefield()
    global minefield minefieldDim
    
    %Could be optimized to only update masked values that have changed
    %(4th page)?
    
    for row = 1:minefieldDim(1)
        for col = 1:minefieldDim(2)
            %If masked, display unknown "-1"
            if(minefield(row, col, 3) && minefield(row, col, 2) ~= 99)
                minefield(row, col, 2) = -1;
            elseif (minefield(row, col, 2) ~= 99)
                %Transfer value from complete minefield
                if(minefield(row, col, 2) ~= minefield(row, col, 1))
                    minefield(row, col, 2) = minefield(row, col, 1);
                    dispCell(row, col); 
                end  
            else
                dispCell(row, col);%display flags
            end
        end
    end
end
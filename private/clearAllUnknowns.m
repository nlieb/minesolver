function clearAllUnknowns()
    global minefield minefieldDim
    
    %Unmask all unknowns
    for i = 1:minefieldDim(1)
        for j = 1:minefieldDim(2)
            if minefield(i,j,2) == -1
                minefield(i,j,3) = 0;
            end
        end
    end
    
    %Apply changes to the masked minefield
    updateMaskedMinefield();
end
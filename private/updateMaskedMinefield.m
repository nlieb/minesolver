%Updates the masked minefield with the latest masking
function updateMaskedMinefield()
    global minefield minefieldDim solvedArray grain    
    
    for m = 1:minefieldDim(1)
        n = 1;
        while n <= minefieldDim(2)
            if(mod(n-1, grain) == 0 && solvedArray(ceil(m/grain),ceil(n/grain)))%dnc optimization
                n = n+grain;%skip to next block
                continue;
            end
            
            %If masked, display unknown "-1"
            if(minefield(m, n, 3) && minefield(m, n, 2) ~= 99)
                minefield(m, n, 2) = -1;
            elseif (minefield(m, n, 2) ~= 99)
                %Transfer value from complete minefield
                if(minefield(m, n, 2) ~= minefield(m, n, 1))
                    minefield(m, n, 2) = minefield(m, n, 1);
                    dispCell(m, n);
                end  
            else
                dispCell(m, n);%display flags
            end
            
            n = n+1;
        end
    end
end
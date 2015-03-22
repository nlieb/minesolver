%divide and conquor methods
function dncRecheckMinefield()
    global minefieldDim solvedArray grain;
    return;
    for m = 1:ceil(minefieldDim(1)/grain)
        for n = 1:ceil(minefieldDim(2)/grain)
            if(~solvedArray(m, n)) %if not solved, check if it's solved now
                solvedArray(m, n) = dncRecheckBlock((m-1)*grain+1,(n-1)*grain+1);
            end
        end
    end
    %fprintf('%d', solvedArray');
    %fprintf('\n');
end

function bSolved = dncRecheckBlock(row,col)
    global minefield minefieldDim grain;
    
    for m = (row-1):(row+grain)
        for n = (col-1):(col+grain)
            if(m > 0 && n > 0 && m <= minefieldDim(1) && n <= minefieldDim(2) ...
                    && minefield(m,n,2) == -1)
                bSolved = 0;
                return;
            end
        end
    end
    
    fprintf('Found solved block at (%d, %d)\n', row, col);
    bSolved = 1;
end


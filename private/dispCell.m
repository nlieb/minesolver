% Update only one cell in the figure window.
function dispCell(i, j)
    global minefield minefieldDim sprites;
    
    % Here be dragons
    if minefield(i, j, 2) == -1
        image([(j-1)*35,j*35],[(minefieldDim(1)-i+1)*35,(minefieldDim(1)-i)*35], sprites{10});
    elseif minefield(i, j, 2) == 99
        image([(j-1)*35,j*35],[(minefieldDim(1)-i+1)*35,(minefieldDim(1)-i)*35], sprites{12});
    elseif minefield(i, j, 2) == 0
        image([(j-1)*35,j*35],[(minefieldDim(1)-i+1)*35,(minefieldDim(1)-i)*35], sprites{11});
    else
        image([(j-1)*35,j*35],[(minefieldDim(1)-i+1)*35,(minefieldDim(1)-i)*35], sprites{minefield(i, j, 2)});
    end
end
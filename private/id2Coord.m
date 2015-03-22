%Converts a position ID to a coordinate on the minefield
function coord = id2Coord(id)
    global minefieldDim
    coord(1) = ceil((id)/minefieldDim(2));
    
    if mod(id,minefieldDim(2)) ~= 0
        coord(2) = mod(id,minefieldDim(2));
    else
        coord(2) = minefieldDim(2);
    end
end
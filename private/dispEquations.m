%Displays all information-containing equations in the solver matrix
function dispEquations()
    global solverMatrix sizeS
    
    for m = 1:sizeS(1)
        flag = false;
        for n = 1:sizeS(2)-1
            if(solverMatrix(m,n) ~= 0)
                coord = id2Coord(n);
                fprintf('[%2d %2d] %d ', coord(1), coord(2), solverMatrix(m,n));
                flag = true;
            end
        end
        if flag == true
            fprintf('%5d\n', solverMatrix(m,sizeS(2)));
        end
    end
    fprintf('\n\n');
end
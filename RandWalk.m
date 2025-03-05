%Homework 6
%Advaita Nair
%UID: 206250044

function [x, y] = RandWalk(x0, y0, BC)
    
    p = rand;
    
    %move up
    if p <= 0.2
        x = x0; y = y0 +1;
        %check for boundary
        if y >= BC(1)
            y = BC(1);
        end
    %move down
    elseif 0.2 < p && p <= 0.4
        x = x0; y = y0 - 1;
        %boundary check
        if y <= BC(2)
            y = BC(2);
        end
    %move left
    elseif 0.4 < p && p <= 0.6
        x = x0 - 1; y = y0;
        %boundary check
        if x <= BC(3)
            x = BC(3);
        end
    %move right
    elseif 0.6 < p && p <= 0.8
        x = x0 + 1; y = y0;
        %boundary check
        if x >= BC(4)
            x = BC(4);
        end
    %no movement
    elseif 0.8 < p
        x = x0; y = y0;
    end
end
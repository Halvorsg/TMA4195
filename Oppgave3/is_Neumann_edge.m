function res = is_Neumann_edge(p1,p2)
    if  (p1(2) == 1 && p2(2) == 1)
        res = 1;
%         fprintf('Vertical\n')
    elseif (p1(2) == -1 && p2(2) == -1)
        res = 2;
%         fprintf('Horizontal\n')
    else
        res = 0;
    end
end
        
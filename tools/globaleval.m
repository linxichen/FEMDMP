function val = globaleval(x,y,xbounds,ybounds,allnodesval)
% Given value of x, y, boundaries of elements, and values at all nodes, use
% local interpolation to find values at (x,y).
% Assuming x,y are within xbounds and ybounds, otherwise return error
if x > xbounds(end)
    error('x too big');
end
if x < xbounds(1)
    error('x too small');
end
if y > ybounds(end)
    error('y too big');
end
if y < ybounds(1)
    error('x too big');
end
i_right = find(xbounds>=x,1,'first');
i_left = max(1,i_right-1);
i_up = find(ybounds>=y,1,'first');
i_down = max(1,i_up-1);
xxi = (x - xbounds(i_left))/max(1,(xbounds(i_right)-xbounds(i_left)));
eeta = (y - ybounds(i_down))/max(1,(ybounds(i_up)-ybounds(i_down)));
nodesval(1) = allnodesval(i_left,i_down);
nodesval(2) = allnodesval(i_right,i_down);
nodesval(3) = allnodesval(i_right,i_up);
nodesval(4) = allnodesval(i_left,i_up);
val = localeval(xxi,eeta,nodesval);
end
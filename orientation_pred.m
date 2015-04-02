function [output] = orientation_pred(ax, ay, bx, by, cx, cy)
% if left then output = 1
% if right then output = 0
d = (ax-cx).*(by-cy) - (bx-cx).*(ay-cy);

output = (d>0);

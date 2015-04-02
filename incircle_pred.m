function [output] = incircle_pred(ax, ay, bx, by, cx, cy, dx, dy)
% if in circle then output = 1
% if not in circle then output = 0
% Note: the order of a,b,c matters. Should be right hand ruled.

d = (ax-dx).*( (by-dy).*((cx-dx).^2+(cy-dy).^2) - (cy-dy).*((bx-dx).^2+(by-dy).^2))...
  - (ay-dy).*( (bx-dx).*((cx-dx).^2+(cy-dy).^2) - (cx-dx).*((bx-dx).^2+(by-dy).^2))...
  + ((ax-dx).^2+(ay-dy).^2).*( (bx-dx).*(cy-dy) - (cx-dx).*(by-dy) );

output = (d>-1e-10);

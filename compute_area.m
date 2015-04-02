function [A] = compute_area(ax,ay,bx,by,cx,cy)

% Compute the area of a triangle
% V1 x V2
% V1 = b - a
% V2 = c - a
%         | bx-ax  by-ay |
% A = 1/2 | cx-ax  cy-ay |

A = 1/2*( (bx-ax).*(cy-ay) - (cx-ax).*(by-ay) );
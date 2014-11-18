function [xi_tri] = transform_quad_to_triangle(xi_quad, range_choice)
% This function contains a mapping a function which trnsforms a 
% quadrilateral into a triangle for use with any quadrature
% Inputs:
%        xi_quad : x and y locations in a quadrilateral ranging from [0,1]
% Outputs:
%        xi_tri : x and y locations in a triangle with vertices [0,0], 
%                 [1,0], 0,1]

% if quad range is [-1,1] 
% from "Quadrature Formulas in Two Dimensions" Math5172 - Finite Element Methods, author?

if (range_choice == 1)
  a = xi_quad(1);
  b = xi_quad(2);
  xi_tri(1) = (1+a).*(1-b)/4;
  xi_tri(2) = (1+b)/2;

elseif (range_choice == 2)
% Same equation except triangle is transformed to a range of [0,1] by substituting a = 2*a_bar - 1 and b = 2*b_bar - 1
% xi_tri(1) = (1+(2a_bar -1))*(1-(2*b_bar-1))/4
% xi_tri(1) = (2*a_bar)*(2-2*b_bar)/4
% xi_tri(1) = 4*(a_bar)*(1-b_bar)/4
% xi_tri(1) = a_bar*(1-b_bar)
%
% xi_tri(2) = (1+(2*b_bar-1))/2
% xi_tri(2) = (2*b_bar)/2
% xi_tri(2) = b_bar

  a_bar = xi_quad(1);
  b_bar = xi_quad(2);
  xi_tri(1) = a_bar.*(1-b_bar);
  xi_tri(2) = b_bar;

end

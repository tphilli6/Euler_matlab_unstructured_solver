function [ index ] = NqD( q, D )
%NqD Calculates all sets of indices in N_d^D in Smolyak quadrature
%   Given q and D, this function calculates all D-tuples such that
%   \sum_{d=1}^D i_d = D + q   (where i_d is the dth index of the D-tuple)
%
%   For simplicity, this routine builds these index sets using recursion.
%%

  index = [];
  count = 0;

  if ( D==2 )  % base case
    for i=1:q+D-1
      index(i,:) = [ i D+q-i ];
    end
  
  else
    for d=1:q+1
      isub = NqD( q+1-d, D-1);
      ni   = size(isub,1);
      for i=1:ni
        count = count + 1;
        index(count,:) = [d isub(i,:)];
      end
    end
  end


end

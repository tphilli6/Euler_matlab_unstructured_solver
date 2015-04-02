function ti = remove_corners(ti)

% Remove the corners
cnt = 1;
for i = 1:size(ti,1)
   if any(ti(i,:)==1) || ...
      any(ti(i,:)==2) || ...
      any(ti(i,:)==3) || ...
      any(ti(i,:)==4) 
  
   else
       ticlean(cnt,:) = ti(i,:);
       cnt = cnt + 1;
   end
end
ti = ticlean-4;
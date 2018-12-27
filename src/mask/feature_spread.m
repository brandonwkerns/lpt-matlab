function array_out = feature_spread(array_in, np)

  s = size(array_in);
  array_last = array_in;
  array_out = array_in;


  if (np > 0)
    for nn = 1:np

      % shift left
      array_shift = circshift(array_last, -1, 3);
      array_shift(:,:,s(3)) = -1;
      array_out = max(array_out, array_shift);

      % shift right
      array_shift = circshift(array_last, 1, 3);
      array_shift(:,:,1) = -1;
      array_out = max(array_out, array_shift);


      % shift down
      array_shift = circshift(array_last, -1, 2);
      array_shift(:,s(2),:) = -1;
      array_out = max(array_out, array_shift);

      % shift up
      array_shift = circshift(array_last, 1, 2);
      array_shift(:,1,:) = -1;
      array_out = max(array_out, array_shift);

      array_last = array_out ;

    end
  end

end %function

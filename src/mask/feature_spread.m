function array_out = feature_spread(array_in, np)

  %% Use convolution to expand the mask "array_in" a radius of np points.
  %% For this purpose, it takes a 3-D array with the first entry being time.

  array_out = array_in;
  s = size(array_in);
  
  [circle_array_x, circle_array_y] = meshgrid(-np:1:np, -np:1:np);
  circle_array_dist = sqrt(circle_array_x.^2 + circle_array_y.^2);
  circle_array_mask = (circle_array_dist < (np + 0.1));
  circle_array_mask = circle_array_mask / sum(sum(circle_array_mask));

  %% Loop over the times.
  %% For each time, use the convolution to "spread out" the effect of each time's field.
  for tt = 1:s(1)
    array_2d = squeeze(array_in(tt,:,:));
    array_2d_new = array_2d;
    unique_values = unique(array_2d);
    unique_values(unique_values == 0) = []; %take out zero -- it is not a feature.
    for this_value = [unique_values]
      starting_mask = (array_2d == this_value);
      starting_mask_spread = conv2(starting_mask, circle_array_mask, 'same');
      array_2d_new(starting_mask_spread > 0.001) = this_value;
    end
    array_out(tt,:,:) = array_2d_new;
  end

end %function

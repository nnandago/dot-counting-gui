function [shift, shifted_im2] = calculate_shift(im1, im2)
%%% calculate how much im2 has shifted relative to im1

    r = size(im2, 1)/4; c = size(im2, 2)/4; height = size(im2, 1)/2; width = size(im2, 2)/2;
    sub_im2 = imcrop(im2, [r c height width]);
    cc = normxcorr2(sub_im2, im1);
    [~, ind_max] = max(abs(cc(:)));
    [max_row, max_col] = ind2sub(size(cc), ind_max);
    
    row_shift = max_row - size(sub_im2, 1) - r; 
    col_shift = max_col - size(sub_im2, 2) - c;
    shift = [row_shift, col_shift];
    
    shifted_im2 = 0*im2;
    row_start = max(1, row_shift); row_end = min(size(im2, 1), size(im2, 1) + row_shift); 
    col_start = max(1, col_shift); col_end = min(size(im2, 2), size(im2, 2) + col_shift);
    shifted_im2(row_start:row_end, col_start:col_end) = im2(1:row_end - row_start + 1, 1:col_end - col_start + 1);
   
end
function dot_thresholds = detect_dot_thresh(imdata, num_channels, seg_im)
   
    nonseg_im = uint16(1 - logical(seg_im));
    dot_thresholds = zeros(1, num_channels); 

    for k = 1:num_channels
        im_stack = imdata{k};
        for p = 1:size(im_stack, 3)
            log_filt_imstack(:, :, p) = medfilt2(im_stack(:, :, p), [3 3]); %logMask(im_stack(:, :, p)).*uint16(logical(seg_im));
        end
        %log_filt_imstack = imgaussfilt3(im_stack); %0*im_stack;
        pot_dots_im = imregionalmax(log_filt_imstack);
        
        max_imstack = max(log_filt_imstack, [], 3);
        frac_thres = graythresh(mat2gray(nonzeros(max_imstack.*nonseg_im)));
        dot_thresholds(k) = frac_thres*max(nonzeros(max_imstack.*nonseg_im));
    end  
end
%         thresh_stack = bsxfun(@times, log_filt_imstack, uint16(seg_im)); thresh_stack(thresh_stack < thresh) = 0;
%         conn_list = bwconncomp(thresh_stack);
%         bright_dot_intensities = [];
%         for p = 1:conn_list.NumObjects
%             if length(conn_list.PixelIdxList{p}) >= 5
%                 bright_dot_intensities = cat(1, bright_dot_intensities, median(log_filt_imstack(ind2sub([512, 512, size(im_stack, 3)], conn_list.PixelIdxList{p}))));
%             else
%                 thresh_stack(ind2sub([512, 512, size(im_stack, 3)], conn_list.PixelIdxList{p})) = 0;
%             end
%         end
%         
%         max_threshstack = max(thresh_stack, [], 3);
%         dilated_im = imdilate(max_threshstack, strel('square', 3));
%         dot_mask = imdilate(dilated_im, strel('square', 3)) - dilated_im;
%         
%         
%         %                 max_imstack(max_imstack < thresh) = 0;
%         %                 frac_thres = graythresh(mat2gray(nonzeros(max_imstack.*nonseg_im)));
%         %                 thresh = frac_thres*max(nonzeros(max_imstack.*nonseg_im));
%         
%         
%         thresh = graythresh(mat2gray(log_filt_imstack(pot_dots_im)))*max(log_filt_imstack(pot_dots_im));
%         good_dots_im = pot_dots_im & log_filt_imstack > thresh;
%         all_dots = logical(sum(good_dots_im, 3)); all_dots = bwlabel(all_dots);
%         if isempty(dot_ims)
%             dot_ims = all_dots;
%         else
%             dot_ims(:, :, end + 1) = all_dots;
%         end
%     end
% 
%     common_dots = logical(dot_ims(:, :, 1).*dot_ims(:, :, 2));
%     for k = 1:size(dot_ims, 3)
%         dots_to_remove = unique(nonzeros(common_dots.*dot_ims(:, :, k)));
%         temp_im = dot_ims(:, :, k);
%         if ~isempty(dots_to_remove)
%             temp_im(temp_im == dots_to_remove) = 0;
%             dot_ims(:, :, k) = bwlabel(temp_im);
%         end
%     end
% end

function lapFrame = logMask(im)   %**check that have the right logMask function

blur_size = 10; 
g_filt = fspecial('gaussian', blur_size, 0.3*(blur_size - 1) + 0.8);

k = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];

lapFrame = imfilter(im, g_filt, 'replicate', 'conv'); %imfilter(im,k,'repl');
end
function dots = detect_dots(imdata, seg_im, num_channels, thresholds)
    
    dots = struct; DOT_AREA = 5; centroids = cell(num_channels, 1);
    for k = 1:num_channels
        im_stack = imdata{k};
        for p = 1:size(im_stack, 3)
            %log_filt_imstack(:, :, p) = medfilt2(im_stack(:, :, p), [3 3]); %logMask(im_stack(:, :, p)).*uint16(logical(seg_im));
        end
        log_filt_imstack = imgaussfilt3(im_stack, 0.8); %0*im_stack;
        pot_dots_im = imregionalmax(log_filt_imstack);
                
        thresh_stack = bsxfun(@times, log_filt_imstack, uint16(seg_im)); thresh_stack(thresh_stack < thresholds(k)) = 0;
        conn_list = bwconncomp(thresh_stack, 26);        
        stats = regionprops(conn_list);
        idx_to_keep = [stats.Area] >= DOT_AREA;
                
        good_dots_stats = stats(idx_to_keep);
        dots(k).properties = good_dots_stats;
        dots(k).counts = length(good_dots_stats);  
        centroids{k} = reshape([good_dots_stats.Centroid], 3, length(good_dots_stats))';
    end
    
    % compare other channels to first channel and eliminate dots
    % that colocalize with dots in the first channel
    OVERLAP_SIZE = 5; dots1_to_remove = [];
    for k = 2:num_channels
        [pairwise_distances, indices] = pdist2(centroids{1}, centroids{k}, 'euclidean', 'Smallest', 1);
        dots_to_keep = find(~ismember(1:length(dots(k).properties), find(pairwise_distances < OVERLAP_SIZE)));
        dots(k).properties = dots(k).properties(dots_to_keep); dots(k).counts = length(dots_to_keep);
        dots1_to_remove = [unique(indices(pairwise_distances < OVERLAP_SIZE))'; dots1_to_remove(:)];
    end
    
    dots1_to_keep = find(~ismember(1:length(dots(1).properties), dots1_to_remove));
    dots(1).properties = dots(1).properties(dots1_to_keep); dots(1).counts = length(dots1_to_keep);
end


% --------------------------------------------------------------------
function lapFrame = logMask(im)   %**check that have the right logMask function

    k = [-4 -1  0 -1 -4;...
         -1  2  3  2 -1;...
          0  3  4  3  0;...
         -1  2  3  2 -1;...
         -4 -1  0 -1 -4];

    lapFrame = imfilter(im,k,'repl');
end
function [seg_im, nuc_im] = segment_on_bg(seg_channel_data, SEG_THRES, nuc_channel_data, NUC_THRES)
%%% segments cells based on autofluorescence in seg_channel. 
%%% if nuclei are fluorescent in separate channel nuc_channel, segments
%%% them too

    %SEG_THRES = 1500; NUC_THRES = 1000;
    num_slices = size(seg_channel_data, 3);

    sum_seg_channel = 0*seg_channel_data(:, :, 1);
    sum_nuc_channel = 0*seg_channel_data(:, :, 1);
    for k = 1:num_slices
        sum_seg_channel = sum_seg_channel + seg_channel_data(:, :, k);
        sum_nuc_channel = sum_nuc_channel + nuc_channel_data(:, :, k);
    end
    
    temp_seg_im = 0*sum_seg_channel; nuc_im = temp_seg_im;
    temp_seg_im(sum_seg_channel > SEG_THRES) = 1;
    nuc_im(sum_nuc_channel > NUC_THRES) = 1;

    se = strel('square', 3);
    temp_seg_im = imopen(temp_seg_im, se); 
    temp_seg_im = imopen(temp_seg_im, se);
    temp_seg_im = imclose(temp_seg_im, se);
    temp_seg_im = imfill(temp_seg_im, 'holes');
    temp_seg_im = bwlabel(temp_seg_im);
    
    nuc_im = imopen(nuc_im, se);
    nuc_im = imopen(nuc_im, se);
    nuc_im = bwlabel(nuc_im, 4);
    props = regionprops(nuc_im);
    areas = [props.Area];
    for k = 1:length(areas)
        if areas(k) < 3000
            nuc_im(nuc_im == k) = 0;
        end
    end
    nuc_im = logical(nuc_im);
    
    cell_nums = unique(nonzeros(temp_seg_im.*nuc_im));
    seg_im = 0*temp_seg_im;
    for k = 1:length(cell_nums)
        seg_im(temp_seg_im == cell_nums(k)) = 1;
    end
    seg_im = bwlabel(seg_im);
end
function [imdata, num_channels] = czi_open(filename)

    im_series = bfopen(filename); 
    num_channels = im_series{4}.getChannelCount(0);
    imdata = cell(num_channels, 1);
    for k = 1:num_channels
        imdata{k} = cat(3, im_series{1}{k:num_channels:end, 1});
    end
end
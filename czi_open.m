function imdata = czi_open(filename, num_channels)

    im_series = bfopen(filename); imdata = cell(num_channels, 1);
    for k = 1:num_channels
        imdata{k} = cat(3, im_series{1}{k:num_channels:end, 1});
    end
end
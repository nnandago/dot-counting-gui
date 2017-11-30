function [ dots, cell_area ] = vectorize_dot_area( stack )
%[ dots, cell_area ] = vectorize_dot_area( stack )
%   Detailed explanation goes here
dots = [];
cell_area = [];
for fov=1:stack.number_images
    if ~isempty(stack.frame(fov).cell)
        for cellNo=1:stack.frame(fov).number_cells
            dots = cat(1,dots,[stack.frame(fov).cell(cellNo).dots.counts]);
            cell_area = cat(1,cell_area,sum(stack.frame(fov).cell(cellNo).mask(:)));
        end
    end
end
end


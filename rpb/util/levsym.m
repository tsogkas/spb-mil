function sym = levsym(img,parts, part_labels, selected_parts_indicator)
% function sym = levsym(img,parts, part_labels, selected_parts_indicator)
% 
% Return binary image with medial axes, as they were detected by
% Levinshtein method for multiscale symmetric parts.

    num_parts = size(parts, 3);
    num_clusters = max(part_labels);
    
    if (nargin < 4 || isempty(selected_parts_indicator))
        selected_parts_indicator = true(num_parts, 1);
    end
                    
    % Draw the medial axes next
    sym = zeros(size(img,1),size(img,2));
    for i = 1:num_parts
        if (selected_parts_indicator(i))
            axis_line = PartMedialAxisLinear(parts(:,:,i));
            x1 = axis_line.y1; x2 = axis_line.y2;
            y1 = axis_line.x1; y2 = axis_line.x2;
            sym = add_segment(sym,x1,y1,x2,y2);
        end
    end
function [ vertices_values ] = faces2vertices( faces_values )
filter = [1 1 ; 1 1]/4;
vertices_values = conv2(faces_values, filter);
end


function [ logicalImg ] = generate_circle_image(logicalImg, R, X, Y )
% creates 2D "logical" Matrix of a circle of radius R, located at (X,Y)
% and adding it to the given logial matrix

[Rows,Cols] = size(logicalImg);
[columnsInImage, rowsInImage] = meshgrid(1:Cols, 1:Rows);

circle = (rowsInImage - Y).^2 ...
    + (columnsInImage - X).^2 <= R.^2;

logicalImg = logicalImg | circle;
end


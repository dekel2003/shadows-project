
%% create new image

% setting image size to 640x480
imageSizeX = 640;
imageSizeY = 480;

logicalImg = false([imageSizeY imageSizeX]);

%% generate 2 circle images I0 and I1
center1 = [30,30]';
center2 = [450,610]';
I0 = generate_circle_image(logicalImg, 20, center1(1), center1(2) );
I1 = generate_circle_image(logicalImg, 20, center2(1), center2(2) );

I0 = I0(:);
I1 = I1(:);

%% required solution
number_of_vertices = numel(I0);
circles_gap = center2 - center1;
v = repmat(circles_gap, number_of_vertices, 1);

%%

laplacian = [];






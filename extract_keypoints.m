function [x, y, scores, Ih, Iv] = extract_keypoints(image)
%% extract_keypoints
% input:
%   image: colored image. Not grayscale or double yet.
% output:
%   x: n x 1 vector of x (col) locations that survive non-maximum suppression. 
%   y: n x 1 vector of y (row) locations that survive non-maximum suppression.
%   scores: n x 1 vector of R scores of the keypoints correponding to (x,y).
%   Ih: x- (horizontal) gradient. Also appeared as Ix in the slides.
%   Iv: y- (vertical) gradient. Also appeared as Iy in the slides.

% The kernels are provided, but you can try other kernels.
Ih_kernel = [1 0 -1; ...
             2 0 -2; ...
             1 0 -1];
Iv_kernel = Ih_kernel';

%% [10 pts] Part A: Setup (Implement yourself)
k = 0.05;
window = 5;
grayIm = double(rgb2gray(image));
Ih = imfilter(grayIm, Ih_kernel);
Iv = imfilter(grayIm, Iv_kernel);
R = zeros(size(grayIm));

%% [15 pts] Part B: R score matrix (Implement yourself)
[rows, cols] = size(grayIm);
for r = 2:rows-1
    for c = 2:cols-1
        M = zeros(2,2);
        for i = r-1:r+1 
            for j  = c-1:c+1
                M(1,1) = M(1,1) + Ih(i,j)^2;
                M(1,2) = M(1,2) + Ih(i,j) * Iv(i,j);
                M(2,1) = M(2,1) + Ih(i,j) * Iv(i,j);
                M(2,2) = M(2,2) + Iv(i,j)^2;
                R(r,c) = det(M) - k * trace(M)^2;
            end
        end
    end
end
R(R==0) = -Inf;


%% Part C: Thresholding R scores (Provided to you, do not modify)
% Threshold standards is arbitrary, but for this assignment, I set the 
% value of the 1th percentile R as the threshold. So we only keep the
% largest 1% of the R scores (that are not -Inf) and their locations.
R_non_inf = R(~isinf(R));
top_R = sort(R_non_inf(:), 'descend');
R_threshold = top_R(round(length(top_R)*0.01));
R(R < R_threshold) = -Inf;

%% [15 pts] Part D: Non-maximum Suppression (Implement yourself)
checkM = zeros(3,3);
index = ones(3,3);
R(1,:) = -Inf;
R(size(R,1),:) = -Inf;
[rows, cols] = size(grayIm);
for i = 2:rows-1
    R(i,1) = -Inf;
    for j  = 2:cols-1
        if (j == cols-1)
            R(i,cols) = -Inf;
        end
        checkM = R(i-1:i+1, j-1:j+1);
        index = (checkM < R(i,j));
        index(2,2) = 1;
        if(any(index(:) == 0))
            R(i,j) = -Inf;
        end
    end
end

x = [];
y = [];
[rows, cols] = size(grayIm);
for i = 1:rows
    for j = 1:cols
        if (R(i,j) ~= -Inf )
            x = [x, j]; % j = col
            y = [y, i]; % i = row
        end
    end
end
scores = zeros(size(x));
for i = 1:length(x)
    scores(i) = R(y(i), x(i));
end
scores = scores';
x = x';
y = y';
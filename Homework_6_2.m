%%
% Load images for Laplacian image blending
img1 = imread('ngc6543a.jpg');
%img2 = imread('street1.jpg');
img2 = imread('rubiks.png');
% Converting images to same size of 600 x 480 pixels
img1 = imresize(img1, [600 480],"bicubic","Antialiasing",true);
img2 = imresize(img2, [600 480],"bicubic","Antialiasing",true);
% Converting images to double for processing
img1 = im2double(img1); 
img2 = im2double(img2); 
% Creating a binary mask from image1 using thresholding
img1_gray = rgb2gray(img1);
mask = imbinarize(img1_gray);
mask = im2double(mask);
% Applying Gaussian filter on the mask for smooth blending
mask = imgaussfilt(mask, 5); % higher sigma for smoothness
% Displaying resized images
fig0 = figure;
figure(fig0);
subplot(3, 1, 1); imshow(img1); title("Nebula - Foreground");
subplot(3, 1, 2); imshow(mask); title("Nebula Mask");
subplot(3, 1, 3); imshow(img2); title("Rubiks - Background");
saveas(fig0, "Images and Mask.png");

% Creating multiresolution pyramid
mrp1 = multiresolutionpyramid(img1, 4);
mrp2 = multiresolutionpyramid(img2, 4);

% Function to create multiresolution pyramid
function mrp= multiresolutionpyramid(A, num_levels)
%   adapted from  https://blogs.mathworks.com/steve/2019/04/09/multiresolution-image-pyramids-and-impyramid-part-2
    M = size(A, 1);
    N = size(A, 2);
    if nargin < 2
        lower_limit = 32;
        num_levels = min(floor(log2([M N]) - log2(lower_limit))) + 1;
    else
        num_levels = min(num_levels, min(floor(log2([M N]))) + 2);
    end
    % Initialize the pyramid structure
    mrp = cell(1, num_levels);
    smallest_size = [M N] / 2^(num_levels - 1);
    smallest_size = ceil(smallest_size);
    padded_size = smallest_size * 2^(num_levels - 1);
    Ap = padarray(A,padded_size - [M N],'replicate','post');
    mrp{1} = Ap;
    for k = 2:num_levels
        mrp{k} = imresize(mrp{k-1},0.5,"bicubic");
    end
    mrp{1} = A;
end

% Visualize the multi resolution pyramid
fig1 = figure;
figure(fig1);
visualizePyramid(mrp1); title("Multi Resolution Pyramid - Nebula - Foreground");
saveas(fig1, "Multi Resolution Pyramid - Nebula - Foreground.png");
fig2 = figure;
figure(fig2);
visualizePyramid(mrp2); title("Multi Resolution Pyramid - Rubiks - Background");
saveas(fig2, "Multi Resolution Pyramid - Rubiks - Background.png");

% Function to visualize the multi resolution pyramid
function tiles_out = visualizePyramid(p)
    % Steve Eddins
    % MathWorks
    M = size(p{1},1);
    N = size(p{1},2);
    for k = 1:numel(p)
        Mk = size(p{k},1);
        Nk = size(p{k},2);
        Mpad1 = ceil((M - Mk)/2);
        Mpad2 = M - Mk - Mpad1;
        Npad1 = ceil((N - Nk)/2);
        Npad2 = N - Nk - Npad1;
        A = p{k};
        A = padarray(A,[Mpad1 Npad1],0.5,'pre');
        A = padarray(A,[Mpad2 Npad2],0.5,'post');
        p{k} = A;
    end
    tiles = imtile(p,'GridSize',[NaN 2],'BorderSize',20,'BackgroundColor',[0.3 0.3 0.3]);
    imshow(tiles)
    if nargout > 0
        tiles_out = tiles;
    end
end

% Applying laplacian pyramid to each image
lapp1 = laplacianPyramid(mrp1);
lapp2 = laplacianPyramid(mrp2);

% Function to generate laplacian pyramid
function lapp = laplacianPyramid(mrp)
    % Steve Eddins
    % MathWorks
 
    lapp = cell(size(mrp));
    num_levels = numel(mrp);
    lapp{num_levels} = mrp{num_levels};
    for k = 1:(num_levels - 1)
    A = mrp{k};
    B = imresize(mrp{k+1},2,"bicubic");
    [M,N,~] = size(A);
    lapp{k} = A - B(1:M,1:N,:);
    end
    lapp{end} = mrp{end};
end

% Visualize laplacian pyramids
fig3 = figure;
figure(fig3);
showLaplacianPyramid(lapp1); title("Laplacian Pyramid - Nebula - Foreground");
saveas(fig3, "Laplacian Pyramid - Nebula - Foreground.png");
fig4 = figure;
figure(fig4);
showLaplacianPyramid(lapp2); title("Laplacian Pyramid - Rubiks - Background");
saveas(fig4, "Laplacian Pyramid - Rubiks - Background.png");

function showLaplacianPyramid(p)
    % Steve Eddins
    M = size(p{1},1);
    N = size(p{1},2);
    stretch_factor = 3;
    for k = 1:numel(p)
        Mk = size(p{k},1);
        Nk = size(p{k},2);
        Mpad1 = ceil((M - Mk)/2);
        Mpad2 = M - Mk - Mpad1;
        Npad1 = ceil((N - Nk)/2);
        Npad2 = N - Nk - Npad1;
        if (k < numel(p))
            pad_value = -0.1/stretch_factor;
        else
            pad_value = 0.4;
        end
        A = p{k};
        A = padarray(A,[Mpad1 Npad1],pad_value,'pre');
        A = padarray(A,[Mpad2 Npad2],pad_value,'post');
        p{k} = A;
    end

    for k = 1:(numel(p)-1)
        p{k} = (stretch_factor*p{k} + 0.5);
    end
    imshow(imtile(p,'GridSize',[NaN 2],'BorderSize',20,'BackgroundColor',[0.3 0.3 0.3]))
end

% Creating a Gaussian pyramid of the mask
gp_mask = gaussianPyramid(mask, 4);

% Function to create gaussian pyramid
function gp = gaussianPyramid(A, num_levels)
% Determining how many iterations are needed
    M = size(A, 1);
    N = size(A, 2);
    if nargin < 2
        lower_limit = 32;
        num_levels = min(floor(log2([M N]) - log2(lower_limit))) + 1;
    else
        num_levels = min(num_levels, min(floor(log2([M N]))) + 2);
    end
    % Initialize the pyramid structure
    gp = cell(1, num_levels);
    gp{1} = A;
    for iter = 2:num_levels
        % Perform gaussian pyramid reduction on the prior level
        gp{iter} = impyramid(gp{iter-1}, 'reduce');
    end
end

% Blending images (Image 1, Mask, Image2)
blended_img = reconstructFromLaplacianPyramid(lapp1, lapp2, gp_mask);

% Function to reconstruct image
function img = reconstructFromLaplacianPyramid(lapp1, lapp2, mask)
    num_levels = numel(lapp2);
    img = repmat(mask{num_levels}, [1, 1, 3]).* lapp1{num_levels} + (1- repmat(mask{num_levels}, [1, 1, 3])).* lapp2{num_levels}; % Blending the top level of pyramid
    for k = num_levels-1:-1:1
        img = imresize(img, 2, 'bicubic');
        blend = repmat(mask{k}, [1, 1, 3]).* lapp1{k} + (1- repmat(mask{k}, [1, 1, 3])).* lapp2{k}; % Blend with the mask
        [M, N, ~] = size(blend);
        img = img(1:M, 1:N, :) + blend;
    end
end

% Display the blended image
figure;
imshow(blended_img);
title('Blended Image');

% Save the blended image to a file
imwrite(blended_img, 'blended_image.png');

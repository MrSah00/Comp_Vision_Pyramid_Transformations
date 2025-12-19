%%
min_corel = 0.5;
 
% load image
I = im2gray(imread('L06 sunflower.png'));
I = imresize(I, 2, "bicubic");
% load template
T = im2gray(imread('L06 sunflower Template.png'));
% Make a smaller template
T = imgaussfilt(T,2);
T = imresize(T, 0.5,'bicubic');
fig1 = figure;
ax1 = axes(fig1);
imshow(I,'Parent', ax1);
fig2 = figure;
figure(fig2);
% While loop till reduced image size is 0.2 times original image size
I_tmp = I;
iter = 0;
while numel(I_tmp) / numel(I) >= 0.2
    iter = iter + 1;
    I_tmp = imgaussfilt(I_tmp,2);
    % calculate and display correlation
    c = normxcorr2(T,I_tmp);
    subplot(4,1,iter); 
    surf(c);
    shading flat;
 
    % get local maxina pixels for correl >0.5
    [yvals,xvals] = find(c > min_corel  & islocalmax(c,1) == 1 & islocalmax(c,2) == 1);
    % display result
    %imshow(I_tmp)
    for x=1:length(yvals)
        fprintf('[%d] [%d] [%d]\n',yvals(x),xvals(x),c(yvals(x),xvals(x)));
        yoffSet = iter*(yvals(x)-size(T,1));
        xoffSet = iter*(xvals(x)-size(T,2));
        rectangle(ax1,'Position', [xoffSet,yoffSet,iter*size(T,2),iter*size(T,1)],'EdgeColor',[1 0 0]); 
    end
    
    I_tmp = imresize(I_tmp, 0.8, "bicubic");
end
% Saving fig1 and fig2
saveas(fig1, 'original_image_with_detections.png');
saveas(fig2, 'correlation_surfaces.png');

% Hough Transform Example.
% JCPatra, Date: 31 Mar 2021
% Ref: Peter Corke Book, p. 184
% A solid square rotated counter-clockwise
% by theta degrees
close all; clear;
%% Read Test Pattern and Rotate

% figure, montage({im, ir}) ; grid on, 
I=(im2gray(imread("Canvas_scripts\295671-transformed.jpeg")));
figure
idisp((I));

%% Edge detection with sobel kernel
vertical_edges = edge(I, 'Sobel', 'vertical');
subplot(1, 2, 1);
imshow(I, []);
title('Original Binary Image');
figure
idisp(vertical_edges);

subplot(1, 2, 2);
imshow(vertical_edges, []);
title('Vertical Edges (Sobel)');


%% nonlocal maxima supression
[H,T,R]=hough(vertical_edges);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(~I,T,R,P,'FillGap',5);

for n=1:numel(lines)
    I=insertShape(I,'line',[lines(n).point1(1),lines(n).point1(2),lines(n).point2(1),lines(n).point2(2)]);
end
figure
idisp(I);
title("Vertical edges")

max_difference=0;
index=0;
for n=1:numel(lines)
    difference=abs(lines(1).point1(1)-lines(n).point1(1));
    if difference>max_difference
        max_difference=difference;
        index=n;
    end
end



point1=lines(1).point1(1);
point2=lines(index).point1(1);
cropped_image=I(:,point2:point1,:);
figure
idisp(cropped_image);
title("Vertically Cropped version");


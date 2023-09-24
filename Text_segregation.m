close all;clear;clc;
%letter_templates=imageDatastore("Template_letters");

%% Imports the image, processes it, removes background, and outputs binary black and white image
I = imread("Car_plate_numbers\295671.jpg");
gs = im2gray(I);
gs = imadjust(gs);
H = fspecial("average",2);
gs = imfilter(gs,H,"replicate");

%Uses morphology to seperate text from rest of bill
SEdisk = strel("disk",7);
Ibg = imclose(gs,SEdisk);
gsSub =  Ibg - gs; 
BW = ~imbinarize(gsSub);

figure;
idisp(BW);

%% Identifies locations of black text by treating letters as blobs 

text=iblobs(BW,'class',1);
parent_numbers = [text.parent];
text_indices = find(parent_numbers ~=0 );
area=text(text_indices).area;
medianarea=median(sort(area));
lines=find(area==medianarea);


% area=zeros(numel(text),1);
% % for n=1:numel(text);
% %     area(n)=text(n).area;
% %     text(n).plot_box('g');
% % end
% 
% for n=1:numel(text)
%     if text(n).parent == 1
%         area(n)=text(n).area;
%         text(n).plot_box('g');
%     end
% end

height=text(lines).vmax - text(lines).vmin;
letter_height=mean(height);


%Test to see whether average letter height obtained from blob analysis works better than original setting
% SEdisk = strel("disk",round(letter_height/2));
SEdisk = strel("disk", round(letter_height/2));
Ibg = imclose(gs,SEdisk);
gsSub =  Ibg - gs; 
BW2 = ~imbinarize(gsSub);

text=iblobs(BW2);

figure;
idisp(BW2);
area=zeros(numel(text),1);
text(text_indices).plot_box('g');

% for n=1:numel(text)
%     if text(n).parent==1
%         text(n)
%         text(n).plot_box('g');
%     end
% end

%% Template matching test using @zncc algorithm to find letter S in bill. IMP: CURRENTLY DOESN'T USE RESULTS FROM iblobs()
% % figure
% % idisp(BW);
% letter=rgb2gray(iread('Template_letters/S.jpg','double'));
% original_height = size(letter, 1);
% original_width = size(letter, 2);
% if mod(original_height,2) ==- 0 || mod(original_width,2) ==0
%     new_height = original_height - 1;
%     new_width = original_width - 1;
% 
%     % Trim the matrix to the new dimensions
%     letter = letter(1:new_height, 1:new_width, :);
% end
% 
% 
% Position=isimilarity(letter,BW,@zncc);
% figure
% idisp(BW);
% [mx,p]=peak2(Position,1,'npeaks',10);
% plot_circle(p,30,'fillcolor','b','alpha',0.3,'edgecolor','none');

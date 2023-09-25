close all;clear;clc;
%letter_templates=imageDatastore("Template_letters");


car_plate="Car_plate_numbers\295671-transformed.jpeg";
letter_text=input("Choose the letter you want\n",'s');
letter_address=strcat("Template_letters\",letter_text,".jpg");


plate=blackandwhite(car_plate);
plate=oddandeven(plate);

edges=icanny(plate);
figure, montage({plate, edges}) ; 
grid on , axis on, title('Car plate, car plate with edges') ;
 
%% Find Hough Transform
h = Hough(edges)
figure, h.show(), grid on ,
lines = h.lines()
idisp(plate) , grid on, title(strcat("Detected ",string(numel(lines))," lines")) ;
h.plot('g')

% nonlocal maxima supression
h1 = Hough(edges, 'suppress', 4)
figure, h.show(), grid on ,
lines = h1.lines()
figure, idisp(plate);
h1.plot('g')

gradient=lines.rho.*tan(deg2rad(lines.theta));
low_slopes=find((-0.1<=gradient) & (gradient<=0.1));
rho_lines=zeros(numel(low_slopes),2);

for n=1:numel(low_slopes)
    rho_lines(n,1)=lines(low_slopes(n)).rho;
    rho_lines(n,2)=low_slopes(n);
end

figure, idisp(plate) ; grid on,

rho_lines=sortrows(rho_lines,1);

disp("Calculating differences")
for n = 1:height(rho_lines)
    v_point1=rho_lines(n,1); 
    v_point2=rho_lines(n+1,1);
    difference(n)=v_point2-v_point1;
    disp(difference(n))
    if n==(height(rho_lines)-1)
        break
    end
end

Nj=find(difference==max(difference));


lines(rho_lines(Nj,2)).plot('r.');
lines(rho_lines(Nj+1,2)).plot('r.');

cropped_image=plate(lines(rho_lines(Nj,2)).rho:lines(rho_lines(Nj+1,2)).rho,:,:);
figure
idisp(cropped_image);
title("Cropped version");

%TODO:Create Hough horizontal edge detection 

plate=morphology(cropped_image,round(lines(rho_lines(Nj+1,2)).rho-lines(rho_lines(Nj,2)).rho));
plate=medfilt2(plate, [4 4]);

figure
idisp(plate);
title("After median filtering");
% for n=1:Nj
%     lines(Nj).plot('r.');
% end

% 
% figure(Name='BW after median filtering');
% idisp(plate);
% title('BW after median filtering')

%% Finding blobs in image, identifying letters, reordering them, 
figure
idisp(plate)
title("Blob imaging")
letters=iblobs(plate);
for n=1:numel(letters)
    if letters(n).touch==0
        letters(n).plot_box('y')
        letter_heights(n)=(letters(n).vmax-letters(n).vmin);
    end
end



letter_height=mode(sort(letter_heights));
%letter_height=weightedAverage(letter_heights);


%plate_letters(n)=[];
position=1;
for n=1:numel(letter_heights)
    if abs(((letter_heights(n)-letter_height)/letter_height)*100) <=1
        plate_letters(position)=letters(n);
        letters(n).plot_box('g');
        position=position+1;
    end
end

plate_letters = sortPlateLettersByValue(plate_letters);
plate_letter_collection=cell(length(plate_letters), 1);

figure
title("Letter")
for n=1:numel(plate_letters)
    plate_letter_collection{n}=plate(plate_letters(n).vmin:plate_letters(n).vmax,plate_letters(n).umin:plate_letters(n).umax,:);
    idisp(plate_letter_collection{n});
    pause(1);
end



%%
letter=blackandwhiteletter(letter_address);
letter=imresize(letter,letter_height/height(letter));
letter=oddandeven(letter);

S = isimilarity(letter, plate, @zncc);
% figure,
% idisp(S, 'colormap', 'jet', 'bar')
% title('zncc between letter and plate');

mi = min(min(S));
mx = max(max(S));
disp('min and max of S');
disp([mi mx]);

% the second option, h = 1 --> w =2h+1, finds max in a 3x3 window
[mx, p] = peak2(S, 1, 'npeaks', 5);
% mx = value of peaks
% p = coordinate (u,v) of corresponding peak.
disp('Five zncc peak values');
disp(mx);

disp('Corresponding u,v coordinates');
disp(p);

figure, idisp(plate);
plot_circle(p, 20, 'edgecolor','b');
plot_point(p, 'sequence', 'bold', 'textsize', 18, 'textcolor', 'r');
title('Five peaks of zncc shown on plate image');


%% Function Definitions

function [image]=blackandwhite(input)
plate=iread(input,'double');
gs = im2gray(plate);
gs = imadjust(gs);
H = fspecial("average",2);
gs = imfilter(gs,H,"replicate");
image=gs;
figure
idisp(plate);
title('before');

%Uses morphology to seperate text from rest of bill
% SEdisk = strel("rectangle",[87 width(plate)]);
% Ibg = imclose(gs,SEdisk);
% gsSub =  Ibg - gs;
% gfilter=kgauss(2);
% gsSub=iconv(gfilter,gsSub);
% BW = ~imbinarize(gsSub);
% BW=medfilt2(BW,[10,10]);
% image=BW;
% figure
% idisp(image);
% title('after');
end

function [image]=morphology(input,height)
gs=input;
SEdisk = strel("rectangle",[height width(input)]);
%Ibg = imopen(gs,SEdisk);
Ibg = imclose(gs,SEdisk);
gsSub =  Ibg - gs;
gfilter=kgauss(2);
gsSub=iconv(gfilter,gsSub);
BW = ~imbinarize(gsSub);
BW=medfilt2(BW,[10,10]);
image=BW;
figure
idisp(image);
title('After morphological operation');
end

function [image]=blackandwhiteletter(input)
plate=iread(input,'double');
gs = im2gray(plate);
gs = imadjust(gs);
H = fspecial("average",2);
gs = imfilter(gs,H,"replicate");
gs=imbinarize(gs);
image=gs;
figure
idisp(plate);
title('before');
end

function [plate_letters] = sortPlateLettersByValue(plate_letters)
    % Extract the values from the objects into a numerical array
    letter_positions = [plate_letters.uc];

    % Sort the values in ascending order and obtain the sorted indices
    [sortedValues, sortedIndices] = sort(letter_positions);

    % Rearrange the objects based on the sorted indices
    plate_letters = plate_letters(sortedIndices);
end


function [your_matrix]=oddandeven(image)
[rows, cols] = size(image);
% Check if either dimension is even
if mod(rows, 2) == 0 || mod(cols, 2) == 0
    % If either dimension is even, make it odd
    if mod(rows, 2) == 0
        rows = rows + 1; % Increase the number of rows by 1
    end
    if mod(cols,2) == 0
        cols = cols + 1; % Increase the number of columns by 1
    end

    % Create a new matrix with the odd dimensions
    new_matrix = ones(rows, cols);

    % Copy the data from the old matrix to the new one
    new_matrix(1:size(image, 1), 1:size(image, 2)) = image;

    % Replace the old matrix in the workspace with the new one
    your_matrix = new_matrix;
else
    your_matrix=image;
end
end

function weightedAvg = weightedAverage(inputArray)
    % Initialize variables to store sum and count for each unique element
    uniqueElements = unique(inputArray);
    sumByElement = zeros(size(uniqueElements));
    countByElement = zeros(size(uniqueElements));
    
    % Calculate sum and count for each unique element
    for i = 1:length(uniqueElements)
        element = uniqueElements(i);
        elementIndices = (inputArray == element);
        sumByElement(i) = sum(elementIndices .* inputArray);
        countByElement(i) = sum(elementIndices);
    end
    
    % Calculate the weighted average
    weightedAvg = sum(sumByElement .* countByElement) / sum(countByElement);
end



% figure(Name='IBlobs BW');
% idisp(BW);
% text=iblobs(BW,'class',1);
% for n = 1:numel(text)
%     text(n).plot_box('g');
% end



% parent_numbers = [text.parent];
% text_indices = find(parent_numbers ~=0 );
% area=text(text_indices).area;
% medianarea=median(sort(area));
% lines=find(area==medianarea);
% 
% 
% height=text(lines).vmax - text(lines).vmin;
% letter_height=mean(height);
% 
% 
% %Test to see whether average letter height obtained from blob analysis works better than original setting
% % SEdisk = strel("disk",round(letter_height/2));
% SEdisk = strel("disk", round(letter_height/2));
% Ibg = imclose(gs,SEdisk);
% gsSub =  Ibg - gs; 
% BW2 = ~imbinarize(gsSub);
% 
% text=iblobs(BW2);
% 
% figure;
% idisp(BW2);
% area=zeros(numel(text),1);
% text(text_indices).plot_box('g');
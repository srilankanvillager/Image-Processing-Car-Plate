close all;clear; clc;
letter_templates=imageDatastore("Template_letters");

text_plate_letters="";
car_plate="Car_plate_numbers\plate7.jpg";

plate=blackandwhite(car_plate);
plate=oddandeven(plate);

original_plate=plate;

edges=icanny(plate);
figure, montage({plate, edges}) ; 
grid on , axis on, title('Car plate, car plate with edges') ;
 
%% Find Hough Horizontal edges
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

horizontalpoint1=lines(rho_lines(Nj,2)).rho;
horizontalpoint2=lines(rho_lines(Nj+1,2)).rho;

%cropped_image=plate(lines(rho_lines(Nj,2)).rho:lines(rho_lines(Nj+1,2)).rho,:,:);
cropped_image=plate(horizontalpoint1:horizontalpoint2,:,:);
figure
idisp(cropped_image);
title("Cropped version");

plate=morphology(cropped_image,round(lines(rho_lines(Nj+1,2)).rho-lines(rho_lines(Nj,2)).rho));
%plate=morphology(cropped_image,300);

plate=medfilt2(plate, [4 4]);

% figure
% idisp(plate);
% title("After median filtering");
% for n=1:Nj
%     lines(Nj).plot('r.');
% end

%% Vertical Hough edge transform

vertical_edges = edge(plate, 'Sobel', 'vertical');

% nonlocal maxima supression
[H,T,R]=hough(vertical_edges);
P  = houghpeaks(H,5,'threshold',ceil(0.6*max(H(:))));
lines = houghlines(vertical_edges,T,R,P,'FillGap',5);

for n=1:numel(lines)
    plate_vertical=insertShape(double(plate),'line',[lines(n).point1(1),lines(n).point1(2),lines(n).point2(1),lines(n).point2(2)],LineWidth=5);
end
figure
idisp(plate_vertical);
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

verticalpoint1=lines(1).point1(1);
verticalpoint2=lines(index).point1(1);

plate=plate(:,min(verticalpoint1,verticalpoint2):max(verticalpoint1,verticalpoint2),:);
figure
idisp(plate);
title("Cropped vertical version");

%% Creating final image for blob analysis
plate=original_plate(horizontalpoint1:horizontalpoint2,min(verticalpoint1,verticalpoint2):max(verticalpoint1,verticalpoint2),:);
plate=morphology(plate,round(abs(horizontalpoint1-horizontalpoint2)));
plate=medfilt2(plate, [4 4]);
figure
idisp(plate);
title("Final image before blob processing");


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


% Choose the number of clusters (K)
K = 4; % You need to decide the number of clusters

% Perform k-means clustering
[idx, centers] = kmeans(letter_heights, K);

% Initialize arrays to store cluster heights and cluster sizes
cluster_heights = zeros(1, K);
cluster_sizes = zeros(1, K);

max_cluster_size = 0; % Initialize maximum cluster size
max_cluster_height = 0; % Initialize maximum cluster height
max_valid_cluster_size = 12;


% Calculate the mean height for each cluster and count cluster sizes
for k = 1:K
    cluster_indices = (idx == k);
    cluster_heights(k) = mean(letter_heights(cluster_indices));
    cluster_sizes(k) = sum(cluster_indices);

    if cluster_sizes(k) < max_valid_cluster_size
        if cluster_sizes(k) > max_cluster_size
            max_cluster_size = cluster_sizes(k);
            max_cluster_height = cluster_heights(k);
        end
    end
end

% Display the estimated letter heights and cluster sizes
fprintf('Estimated Letter Heights and Cluster Sizes:\n');
for k = 1:K
    fprintf('Cluster %d: Height=%.2f, Size=%d\n', k, cluster_heights(k), cluster_sizes(k));
end

%letter_height=cluster_heights(find(cluster_heights==max(cluster_heights(find(cluster_sizes<9)))));
letter_height=max(cluster_heights(find(cluster_sizes<9)));
%letter_height=300;

figure
idisp(plate)

plate_letters=RegionFeature.empty(0);


position=1;
for n=1:numel(letter_heights)
    letter_centre=[letters(n).uc, letters(n).vc];
    if abs(((letter_heights(n)-letter_height)/letter_height)*100) <=50.0
        plate_letters(position)=letters(n);
        letters(n).plot_box('g');
        position=position+1;
    end
end

plate_letters = sortPlateLettersByValue(plate_letters);

for n=1:numel(plate_letters)
    possible_letter=plate_letters(n);
    if n==1
        plate_letters(n)=possible_letter;
    end
    if n~=1
        percentage_difference= abs((possible_letter.uc-plate_letters(n-1).uc)/plate_letters(numel(plate_letters)).uc)*100;
        if percentage_difference<=7
            disp("Same letter")
            if possible_letter.area > plate_letters(n-1).area
                plate_letters(n-1)=[];
            end
            if possible_letter.area < plate_letters(n-1).area
                plate_letters(n)=[];
            end
        end
    end
    if n >= numel(plate_letters)
        break
    end
end

nonEmptyIndices = ~arrayfun(@isempty, plate_letters);

% Resize the array to remove empty elements
plate_letters = plate_letters(nonEmptyIndices);


plate_letter_collection=cell(length(plate_letters), 1);

figure
title("Letter")
for n=1:numel(plate_letters)
    plate_letter_collection{n}=plate(plate_letters(n).vmin:plate_letters(n).vmax,plate_letters(n).umin:plate_letters(n).umax,:);
    idisp(plate_letter_collection{n});
    pause(1);
end


%% Iterating over every letter in the blob image matrix to find template letter with highest match

matches=findBestTemplateMatchWithResize(plate_letter_collection,letter_templates);
text="";
for n=1:numel(matches)
    text=strcat(text,stringsplitter(matches{n}));
end

disp(text);

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



function plate_letter_positions = findBestTemplateMatchWithResize(plate_letter_collection, letter_templates)
plate_letter_positions=cell(numel(plate_letter_collection),1);
    for n=1:numel(plate_letter_collection)
        disp("---")
        currentImage = double(plate_letter_collection{n});
        bestMatchScore = inf;  % Initialize with negative infinity
        bestMatchIndex = -1;
        for m=1:numel(letter_templates.Files)
            letter_template=double(rgb2gray(iread(letter_templates.Files{m})));
            letter_template=imresize(letter_template,size(currentImage));
            similarityScore = immse(currentImage, letter_template);
            idisp(similarityScore)
            if similarityScore < bestMatchScore
                bestMatchScore = similarityScore;
                bestMatchIndex = m;
                disp(letter_templates.Files{m});
            end
        end
        plate_letter_positions{n,1}=letter_templates.Files{bestMatchIndex};
    end
end

function letter=stringsplitter(address)
    % Define the input file path as a string
    filePath = 'C:\Users\PC\Desktop\MATLAB minor project\Template_letters\N.jpg';
    
    % Use fileparts to get the base name of the file
    [~, baseFileName, ~] = fileparts(filePath);
    
    % Split the base file name using the file separator (assuming '\' or '/')
    fileNameParts = strsplit(baseFileName, filesep);
    
    % The last part of the split string is the file name
    letter = fileNameParts{end};
end

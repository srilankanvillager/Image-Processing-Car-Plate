% Hough Transform Example.
% JCPatra, Date: 31 Mar 2021
% Ref: Peter Corke Book, p. 184
% A solid square rotated counter-clockwise
% by theta degrees
%% Read Test Pattern and Rotate
im = testpattern('squares', 256, 256, 128);
% ir = imrotate(im,30,'bilinear'); % MATLAB 
theta = -30 ; % i degrees
ir  = irotate(im, theta );  
figure, montage({im, ir}) ; grid on, 

%% Find edges using Canny Algorithm
edges = icanny(ir);
figure, montage({im,ir, edges}) ; 
grid on , axis on, title('sq, rot sq, edgesCanny') ;
% figure, idisp(1-edges), grid on ; axis on;


%% Find Hough Transform
h = Hough(edges)
figure, h.show(), grid on ,
lines = h.lines()
% shows all 9 detected lines some are very close to each other
idisp(ir) , grid on, title('Detected 9 lines') ;
h.plot('g')

%% nonlocal maxima supression
h1 = Hough(edges, 'suppress', 5)
figure, h.show(), grid on ,
lines = h1.lines()
% Now we get only the FOUR dominant lines.
% The detected lines can be projected onto 
% the original image
figure, idisp(ir);
h1.plot('g')

%% See specific  lines
figure, idisp(ir) ; grid on,
% Display first Nj major lines
Nj = 2 ; % change Nj values and see
lines(1:Nj).plot('r.');


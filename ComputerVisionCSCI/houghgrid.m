grid = imread('grids.png');

grey = rgb2gray(grid);

dubs = im2double(grey);

imbw = im2bw(dubs, 0.8);

BW = 1-imbw;

[H,T,R]  = hough(BW);

P  = houghpeaks(H,25,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);

max_len = 0;
figure(1)
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   hold on
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
   hold on

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end


% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
hold off

figure(3)
imshow(grey)

figure(4)
imshow(BW)


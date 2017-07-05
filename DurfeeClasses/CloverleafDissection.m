Cloverleadf


%%

%Clear all previous windows
clear all;
close all;
clc;

vid = videoinput('winvideo', 1, 'RGB24_320x240');

disp('Beggining Color Tracking')
disp('Press Stop Button To Exit')

camim = getsnapshot(vid);

%Creating the figure window to display video and tracking
f=figure('Name', 'Webcam Capture');
%Adding a stop button to the figure
u=uicontrol('string','stop','callback','delete(get(gcbo,''parent''))');
%Adding a text box to display the area of the object
hText1 = uicontrol('Style','text','Parent',f,...
                   'Position',[20 125 80 80],...
                   'BackgroundColor',[0.7 0.7 0.7]);
%adding a text bow to display distance to object
hText2 = uicontrol('Style','text','Parent',f,...
   'Position',[20 250 80 80],...
   'BackgroundColor',[0.7 0.7 0.7]);
%Displaying the top sub plot will be rgb and the bottom the depth image

subplot(5,5,1:10),h1=imshow(camim); 
subplot(5,5,16:25),h2=imshow(camim,[0 9000]); colormap('jet');
subplot(5,5,13),text(0,0,'Wecam Capture');
set(gca,'Visible','off')


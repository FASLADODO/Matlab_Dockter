%% Plot in 3d, but with patches so that they can be transparent 
%   Format:
% 
%   h = plot3transparent(d, facecolor, alpha, ptsize, cubesInsteadOfTetrahedra)
function h = plot3transparent(d, facecolor, alpha, ptsize, cubesInsteadOfTetrahedra)
% t=[0:.1:2*pi]';
% d= [sin(t) cos(t) t];
% alpha =.1;
% facecolor = 'r'; % [ 1  0 0 ];
% ptsize=.1;
% cubesInsteadOfTetrahedra=1;


%% generate basic 3d object template. 
Tt = [... % tetrahedron
   -1.0000   -0.5000   -0.5000
    1.0000   -0.5000   -0.5000
         0    0.5000   -0.5000
         0         0    0.5000]*ptsize(1,1);
     
Ttf= [3 1 2;  4 2 1 ; 1 3 4; 2 3 4]; % tetrahedron verticies

Cb = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1 ; 1 0 1; 1 1 1 ; 0 1 1];%*ptsize;
%scale each dimension
msz = range(d)./max(range(d));
trPtSize = repmat(msz,8,1);
Cb = (Cb - .5).*(trPtSize*ptsize);
Cbf= [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    
if cubesInsteadOfTetrahedra
    Mrk = Cb;
    Mrkf = Cbf;
else
    Mrk = Tt;%Cb;
    Mrkf = Ttf;%Cbf;
end
msize=size(Mrk,1);
vsize=size(Mrkf,1);

dd = zeros( msize*size(d,1), size(d, 2));
ddf= zeros( vsize*size(d,1), size(Mrkf, 2));

%% generate '3d' markers from data...
for i = 1:length(d)
   for j=1:size(d,2)
      dd((i-1)*msize +[1:msize],j) =  d(i,j) + Mrk(:,j);      
   end
   ddf((i-1)*vsize +[1:vsize],:) = Mrkf + msize*(i-1);
end


%figure
v=size(ddf, 1); % vertex count
Clr =[ones(v,1)*.5 zeros(v,1), zeros(v,1)];
%patch ('Vertices',T','Faces',Tf);
%patch(T(1,:),T(2,:),T(3,:),'r')
h=patch ('Vertices',dd,'Faces',ddf,...
 'FaceColor',facecolor,...
'facealpha',alpha,'edgealpha', 0);
%    'FaceVertexCData','r','FaceColor','flat')

%xlabel('x [cm]'); ylabel('y [cm]'); zlabel('z [cm]')
%grid on
view(3); 
%axis square



return
figure
v=size(Tf, 1); % vertex count
%patch ('Vertices',T','Faces',Tf);
%patch(T(1,:),T(2,:),T(3,:),'r')
patch ('Vertices',T,'Faces',Tf,...
 'FaceVertexCData',[ones(v,1)*.5 zeros(v,1), zeros(v,1)],'FaceColor','flat',...
'facealpha',.5,'edgealpha', .5)
%    'FaceVertexCData','r','FaceColor','flat')

xlabel('x [cm]'); ylabel('y [cm]'); zlabel('z [cm]')
grid on

%%
figure
v=size(Cbf, 1); % vertex count
%patch ('Vertices',T','Faces',Tf);
%patch(T(1,:),T(2,:),T(3,:),'r')
h = patch ('Vertices',Cb,'Faces',Cbf,...
 'FaceVertexCData',[ones(v,1)*.5 zeros(v,1), zeros(v,1)],'FaceColor','flat',...
'facealpha',.5,'edgealpha', .5)
%    'FaceVertexCData','r','FaceColor','flat')

xlabel('x [cm]'); ylabel('y [cm]'); zlabel('z [cm]')
grid on
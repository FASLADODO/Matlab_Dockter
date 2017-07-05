% plot PDF over astrini grasp data
clear all

load segdata.mat

fontS = 14;

state_d = segData;

CM = { [1,0,0] , [0,1,0], [0,0,1], [1,1,0], [0,1,1], [1,0,1], [0,0,0], [1,1,1] };
typearr = [1,1,2,2];
enumTypes = {'gallbladder'; 'smallbowel' };

for kk = 1:max(typearr)
    combine{kk}.angle = [];
    combine{kk}.angledot = [];
    combine{kk}.angledotdot = [];
end

for kk = 1:4
    combine{typearr(kk)}.angle = [combine{typearr(kk)}.angle, state_d{kk}.angle'];
    combine{typearr(kk)}.angledot = [combine{typearr(kk)}.angledot, state_d{kk}.angledot'];
    combine{typearr(kk)}.angledotdot = [combine{typearr(kk)}.angledotdot, state_d{kk}.angledotdot'];
end


%% Plot all the phase and grid and avgs

figure(1)

for kk = 1:max(typearr)
    h(kk) = quiver(combine{kk}.angle,combine{kk}.angledot,combine{kk}.angledot,combine{kk}.angledotdot,'color',CM{kk});

    hold on
end


window = 7;
minx = 1.5;
maxx = 4.5;
miny = 0;
maxy = 11;

tickx1 = (maxx - minx)/window;
tickx2 = (maxy - miny)/window;

x1grid = linspace(minx,maxx,window+1);
x2grid = linspace(miny,maxy,window+1);

% Horizontal grid 
for k = 1:length(x2grid)
  line([x1grid(1) x1grid(end)], [x2grid(k) x2grid(k)])
end
hold on
% Vertical grid
for k = 1:length(x2grid)
  line([x1grid(k) x1grid(k)], [x2grid(1) x2grid(end)])
end
hold on


%plot means
tt = 1;
for kk = 1:max(typearr)
    for ii = 1:length(x1grid)-1
       for jj = 1:length(x2grid)-1
           %segment region
          i1 = find( (combine{kk}.angle >= x1grid(ii) & combine{kk}.angle <= x1grid(ii+1)) & (combine{kk}.angledot >= x2grid(jj) & combine{kk}.angledot <= x2grid(jj+1)) );
          indexer{ii}.ind{jj} = i1;
          
          %find mean in each region
          region{kk}.hor{ii}.vert{jj}.angle = combine{kk}.angle(i1);
          region{kk}.hor{ii}.vert{jj}.avg_angle = mean(combine{kk}.angle(i1));
          region{kk}.hor{ii}.vert{jj}.angledot = combine{kk}.angledot(i1);
          region{kk}.hor{ii}.vert{jj}.avg_angledot = mean(combine{kk}.angledot(i1));
          
          if(~isempty(region{kk}.hor{ii}.vert{jj}.avg_angle) & ~isempty(region{kk}.hor{ii}.vert{jj}.avg_angledot) )
              if(kk == 1)
                plot(region{kk}.hor{ii}.vert{jj}.avg_angle,region{kk}.hor{ii}.vert{jj}.avg_angledot,'bx','markersize',10);
              end
              if(kk == 2)
               plot(region{kk}.hor{ii}.vert{jj}.avg_angle,region{kk}.hor{ii}.vert{jj}.avg_angledot,'kx','markersize',10);
              end
          end
       end
    end
end

axis([min(x1grid),max(x1grid),min(x2grid),max(x2grid)])

% title('Phase Portrait regions (angle, angledot)','FontSize',fontS)
% xlabel('Angle (rad)','FontSize',fontS)
% ylabel('Angledot (rad/s)','FontSize',fontS)
% legend([h(1),h(2)],enumTypes{1},enumTypes{2},'FontSize',fontS)
% 

%% Get all PDFs in each region


count = 1;
lt = 100;
for ii = 1:length(x1grid)-1
    for jj = 1:length(x2grid)-1
%       if(~isempty(region{1}.hor{ii}.vert{jj}.angle) && ~isempty(region{1}.hor{ii}.vert{jj}.angledot) && ~isempty(region{2}.hor{ii}.vert{jj}.angle) && ~isempty(region{2}.hor{ii}.vert{jj}.angledot)  )
        if(length(region{1}.hor{ii}.vert{jj}.angle) > lt && length(region{1}.hor{ii}.vert{jj}.angledot) > lt && length(region{2}.hor{ii}.vert{jj}.angle) > lt && length(region{2}.hor{ii}.vert{jj}.angledot) > lt  )
            %get between class
            data = [region{1}.hor{ii}.vert{jj}.angle,region{2}.hor{ii}.vert{jj}.angle; region{1}.hor{ii}.vert{jj}.angledot,region{2}.hor{ii}.vert{jj}.angledot];
            [bandwidth_bc,density_bc,X_bc,Y_bc]=kde2d(data');
            
            %Get within class
            data_wc1 = [region{1}.hor{ii}.vert{jj}.angle;region{1}.hor{ii}.vert{jj}.angledot];
            data_wc2 = [region{2}.hor{ii}.vert{jj}.angle;region{2}.hor{ii}.vert{jj}.angledot];
            [bandwidth_wc1,density_wc1,X_wc1,Y_wc1]=kde2d(data_wc1');
            [bandwidth_wc2,density_wc2,X_wc2,Y_wc2]=kde2d(data_wc2');
            
%             contour3(X_wc1,Y_wc1,density_wc1,50);
%             hold on
%             contour3(X_wc2,Y_wc2,density_wc2,50);
            
            %store all the things
            region{1}.hor{ii}.vert{jj}.X_bc = X_bc;
            region{1}.hor{ii}.vert{jj}.Y_bc = Y_bc;
            region{1}.hor{ii}.vert{jj}.density_bc = density_bc;
            region{1}.hor{ii}.vert{jj}.X_wc = X_wc1;
            region{1}.hor{ii}.vert{jj}.Y_wc = Y_wc1;
            region{1}.hor{ii}.vert{jj}.density_wc = density_wc1;
            region{2}.hor{ii}.vert{jj}.X_wc = X_wc2;
            region{2}.hor{ii}.vert{jj}.Y_wc = Y_wc2;
            region{2}.hor{ii}.vert{jj}.density_wc = density_wc2;
            %compute entropys
            entropy_bc = -sum(sum( density_bc * log(abs(density_bc+(density_bc==0))))); %ignore elements=0, they screw stuff up
            entropy_wc = -sum(sum( density_wc1* log(abs(density_wc1+(density_wc1==0))))) - sum(sum( density_wc2* log(abs(density_wc2+(density_wc2==0)))));
            %Save weights
            if(~isnan(entropy_bc) || ~isnan(entropy_wc))
                weights(ii,jj) = entropy_wc / entropy_bc;
            else
                weights(ii,jj) = 0;
            end
            count = count + 1;
        else
            region{1}.hor{ii}.vert{jj}.X_bc = [];
            region{1}.hor{ii}.vert{jj}.Y_bc = [];
            region{1}.hor{ii}.vert{jj}.density_bc = [];
            region{1}.hor{ii}.vert{jj}.X_wc = [];
            region{1}.hor{ii}.vert{jj}.Y_wc = [];
            region{1}.hor{ii}.vert{jj}.density_wc = [];
            region{2}.hor{ii}.vert{jj}.X_wc = [];
            region{2}.hor{ii}.vert{jj}.Y_wc = [];
            region{2}.hor{ii}.vert{jj}.density_wc = [];
            weights(ii,jj) = 0;
        end
        
    end
end


% hold off
% title('Phase Portrait + within-class density (angle, angledot)','FontSize',fontS)
% xlabel('Angle (rad)','FontSize',fontS)
% ylabel('Angledot (rad/s)','FontSize',fontS)
% legend([h(1),h(2)],enumTypes{1},enumTypes{2},'FontSize',fontS)

%% Plot entropy ratio weights
%Thanks Internet!
%http://www.mathworks.com/help/matlab/visualize/representing-a-matrix-as-a-surface.html

hold on

cc = 1;
for ii = 1:length(x1grid)-1
    for jj = 1:length(x2grid)-1
        x1plot(cc) = x1grid(ii) + (tickx1/2);
        x2plot(cc) = x2grid(jj) + (tickx2/2);
        he(cc) = weights(ii,jj);
        cc= cc +1;
    end
end

x1lin = linspace(min(x1plot),max(x1plot),33);
x2lin = linspace(min(x2plot),max(x2plot),33);

[X,Y] = meshgrid(x1lin,x2lin);
f = scatteredInterpolant(x1plot',x2plot',he');
Z = f(X,Y);

mesher = mesh(X,Y,Z) %interpolated
colormap cool

hold off
view(0, 20);
set(mesher,'facecolor','none')
title('Phase Portrait with seperability weights','FontSize',fontS)
xlabel('Angle (rad)','FontSize',fontS)
ylabel('Angledot (rad/s)','FontSize',fontS)
zlabel('Entropy Weights','FontSize',fontS)
h_legend=legend([h(1),h(2)],enumTypes{1},enumTypes{2})
set(h_legend,'FontSize',12);

%%

% figure(2)
% 
% for kk = 1:numFiles
%     g(kk) = quiver(segData{kk}.strain,segData{kk}.straindot,segData{kk}.straindot,segData{kk}.straindotdot,'color',CM{FullData{kk}.tissuetype})
% 
%     hold on
% end
% hold off
% 
% title('Phase Portrait (strain, straindot)','FontSize',fontS)
% xlabel('strain (N)','FontSize',fontS)
% ylabel('straindot (N/s)','FontSize',fontS)
% 
% legend([g(typeIndex(1)),g(typeIndex(2))],enumTypes{1},enumTypes{2})
% 
% figure(numFiles + 3)
% 
% for kk = 1:numFiles
%     f(kk) = quiver(segData{kk}.strain,segData{kk}.angledot,segData{kk}.straindot,segData{kk}.angledotdot,'color',CM{FullData{kk}.tissuetype})
% 
%     hold on
% end
% hold off
% 
% title('Phase Portrait (strain, angledot)','FontSize',fontS)
% xlabel('strain (N)','FontSize',fontS)
% ylabel('angledot (rad/s)','FontSize',fontS)
% 
% legend([f(typeIndex(1)),f(typeIndex(2))],enumTypes{1},enumTypes{2})
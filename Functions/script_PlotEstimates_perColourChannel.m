%% This script is used to plot the colour channels separately (Fig. 2c).
% It plots the shape-only estimates for red, green, and blue colour
% channels as separate plots, an example is shown in manuscript Fig. 2c

% Last Modified by John Murray-Bruce at University of South Florida
% 30-Nov-2023 [Clean-up and commented for publishing]

% Define red, green, and blue colour maps
nCols = 100;
map_r = zeros(nCols, 3);
map_g = map_r;
map_b = map_r;
map_r(:,1) = linspace(0,1,nCols);
map_g(:,2) = linspace(0,1,nCols);
map_b(:,3) = linspace(0,1,nCols);

% Plot shape of red channel only (using red colour map)
figure(2+11); clf;
Col_r = sceneShape_col(:,1); 
scatter3(scene_est_xyz(:,1),scene_est_xyz(:,2),scene_est_xyz(:,3),10,Col_r,'filled');
hold on;
axis equal
xlim([0 L_max])
ylim([0 L_max])
zlim([0 L_max])

set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4)
xlabel('$x$','FontSize',font_size,'interpreter','latex')
ylabel('$y$','FontSize',font_size,'interpreter','latex')
zlabel('$z$','FontSize',font_size,'interpreter','latex')
set(gca,'zdir','reverse','TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
% ax = gca; ax.GridLineStyle = '--'; ax.GridColor = 10*[0.1, 0.1, 0.1]; grid on
grid on;
grid minor
box on;
colormap(map_r);

% Plot shape of green channel only (using green colour map)
figure(2+12); clf;
Col_g = sceneShape_col(:,2);
scatter3(scene_est_xyz(:,1),scene_est_xyz(:,2),scene_est_xyz(:,3),10,Col_g,'filled');
hold on;
axis equal
xlim([0 L_max])
ylim([0 L_max])
zlim([0 L_max])

set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4)
xlabel('$x$','FontSize',font_size,'interpreter','latex')
ylabel('$y$','FontSize',font_size,'interpreter','latex')
zlabel('$z$','FontSize',font_size,'interpreter','latex')
set(gca,'zdir','reverse','TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
grid on;
grid minor
box on;
colormap(map_g);


% Plot shape of blue channel only (using blue colour map)
figure(2+13); clf;
Col_b = sceneShape_col(:,3);
scatter3(scene_est_xyz(:,1),scene_est_xyz(:,2),scene_est_xyz(:,3),10,Col_b,'filled');
hold on;
axis equal
xlim([0 L_max])
ylim([0 L_max])
zlim([0 L_max])

set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4)
xlabel('$x$','FontSize',font_size,'interpreter','latex')
ylabel('$y$','FontSize',font_size,'interpreter','latex')
zlabel('$z$','FontSize',font_size,'interpreter','latex')
set(gca,'zdir','reverse','TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
colormap(map_b);
grid on;
grid minor
box on;

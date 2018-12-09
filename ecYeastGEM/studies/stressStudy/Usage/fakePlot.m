%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create data:
levels      = 0:3;
enzyme_up   = (1:1:4) + randn(1,4)/5;
enzyme_down = (5:-5/4:1) + randn(1,4)/5;

%Linear fits:
lmodel_up = fitlm(levels,enzyme_up);
disp(['Fit for lmodel_up: R^2 = ' num2str(lmodel_up.Rsquared.Ordinary)])
lmodel_down = fitlm(levels,enzyme_down);
disp(['Fit for lmodel_down: R^2 = ' num2str(lmodel_down.Rsquared.Ordinary)])

%Plot:
figure
hold on
ms = 10;
lw = 3;
plot(levels,enzyme_up,'o','MarkerEdgeColor','m','MarkerSize',ms,'LineWidth',lw)
plot(levels,enzyme_down,'o','MarkerEdgeColor','c','MarkerSize',ms,'LineWidth',lw)
plot(levels',predict(lmodel_up,levels'),'-m','LineWidth',lw)
plot(levels',predict(lmodel_down,levels'),'-c','LineWidth',lw)
y_lab  = 'Enzyme usage [%]';
x_lab  = 'Stress level';
setOptions(x_lab,[0 3],1:2,y_lab,[0.5 5.5],1:5)
axis square
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

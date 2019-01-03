%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotCorr(x,y,x_lab,x_tick)

%Replace any zero with NaN:
[m,n] = size(x);
for i = 1:n
    if ~contains(x_lab,'stress')
        pos = (x(:,i) == 0) + (y(:,i) == 0) > 0;
    else
        pos = (y(:,i) == 0) > 0;
    end
    x(pos,i) = NaN;
    y(pos,i) = NaN;
end

%Plot unfeasible region:
hold on
if max(max(y)) > 101 && m > 1
    patch([x_tick(1),x_tick(end),x_tick(1)], ...
          [x_tick(1),x_tick(end),x_tick(end)],[0.7 0.7 0.7]);
    text(x_tick(1)*2,x_tick(end-1),'Unfeasible region','FontSize',15)
end

%Plot data:
colors = sampleCVDmap(6);
if n == 1
    colors = 'c';
elseif n == 14
    colors = [[0 0 0]; repmat(colors(6,:),3,1); repmat(colors(2,:),7,1); ...
              repmat(colors(4,:),3,1)];
else
    if contains(x_lab,'Temperature')
        colors = [[0 0 0]; repmat(colors(6,:),3,1)];
    elseif contains(x_lab,'NaCl')
        colors = [[0 0 0]; repmat(colors(2,:),7,1)];
    else
        colors = [[0 0 0]; repmat(colors(4,:),3,1)];
    end
end
x_r = [];
y_r = [];
if m == 1
    ms = 10;
    lw = 2;
else
    ms = 4;
    lw = 0.35;
end
if n == 14
    for i = [5:11 1 2 12 3 13 4 14]
        pos = isnan(x(:,i)+y(:,i)) == 0;
        plot(x(pos,i),y(pos,i),'o','MarkerEdgeColor',colors(i,:),'MarkerSize',ms,...
            'LineWidth',lw)
        x_r = [x_r;x(pos,i)];
        y_r = [y_r;y(pos,i)];
    end
else
    for i = 1:n
        pos = isnan(x(:,i)+y(:,i)) == 0;
        plot(x(pos,i),y(pos,i),'o','MarkerEdgeColor',colors(i,:),'MarkerSize',ms,...
            'LineWidth',lw)
        x_r = [x_r;x(pos,i)];
        y_r = [y_r;y(pos,i)];
    end
end

%Log scales and trendlines:
set(gca,'yscale','log')
if ~contains(x_lab,'MW') && ~contains(x_lab,'stress')
    set(gca,'xscale','log')
    lmodel = fitlm(log10(x_r),log10(y_r));
    plot(x_tick',10.^predict(lmodel,log10(x_tick')),'-m','LineWidth',lw)
else
    lmodel = fitlm(x_r,log10(y_r));
    plot(x_tick',10.^predict(lmodel,x_tick'),'-m','LineWidth',lw)
end
disp(['Fit for ' x_lab ': R^2 = ' num2str(lmodel.Rsquared.Ordinary)])

%Other options:
if max(max(y)) == 100
    y_tick = 10.^(-6:2:2);
    if n == 1
        y_lab  = 'Average usage [%]';
    else
        y_lab  = 'Enzyme usage [%]';
    end
elseif m == 1
    if contains(x_lab,'stress')
        if lmodel.Coefficients{2,1} > 0
            y_lab  = 'Enzyme X usage [%]';
            y_tick = ceil(10.^(1.4:0.05:1.55));
        else
            y_lab  = 'Enzyme Y usage [%]';
            y_tick = ceil(10.^(0.6:0.15:1.05));
        end
    else
        y_lab  = [];
        y_tick = 10.^(-2:3:5);
    end
else
    y_lab  = 'Enzyme usage [nmol/gDW]';
    y_tick = 10.^(-7:2:3);
end
setOptions(x_lab,[x_tick(1) x_tick(end)],x_tick, ...
            y_lab,[y_tick(1) y_tick(end)],y_tick)
axis square
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

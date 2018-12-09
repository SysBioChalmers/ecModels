%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function setOptions(x_lab,x_lim,x_ticks,y_lab,y_lim,y_ticks)

%X-Axis:
text_size = 13;
if ~isempty(x_lab)
    xlabel(x_lab,'FontSize',text_size);
end
if ~isempty(x_lim)
    xlim(x_lim)
end
if ~isempty(x_ticks)
    set(gca,'XTick',x_ticks)
end

%Y-Axis:
if ~isempty(y_lab)
    ylabel(y_lab,'FontSize',text_size);
end
if ~isempty(y_lim)
    ylim(y_lim)
end
if ~isempty(y_ticks)
    set(gca,'YTick',y_ticks)
end

%Other:
set(gca,'FontSize',text_size)
set(gca,'Layer','bottom')
set(gca,'XColor','k')
set(gca,'YColor','k')
box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

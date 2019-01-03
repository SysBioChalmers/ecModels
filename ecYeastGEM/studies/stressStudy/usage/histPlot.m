%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function histPlot(data,N,x_lab,y_lab)

if isempty(N)
    N = ceil(length(data)/20);
end

[counts,centers] = hist(data,N);
if contains(x_lab,'slope')
    bar(centers(counts > 0),counts(counts > 0),'c','BaseValue',0.7)
    set(gca,'yscale','log')
    y_lim   = [0.7 1e3];
    x_lim   = [-10,10];
    y_ticks = 10.^(0:3);
else
    bar(centers,counts,'c','BarWidth',1)
    y_lim   = [];
    x_lim   = [];
    y_ticks = [];
end
if N > 50
    set(get(gca,'child'),'EdgeColor','none');
end
setOptions(x_lab,x_lim,[],['Number of ' y_lab],y_lim,y_ticks)
if contains(x_lab,'slope')
    set(gca,'YTickLabel',{'1','10','100','1000'});
end
axis square

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

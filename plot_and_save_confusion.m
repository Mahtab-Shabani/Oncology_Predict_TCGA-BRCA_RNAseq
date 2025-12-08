function plot_and_save_confusion(mat_counts, order, titleStr, fname)
    fig = figure('Visible','off','Position',[200 200 560 420]);
    imagesc(mat_counts);
    colormap('parula'); colorbar;
    axis equal tight;
    set(gca,'XTick',1:length(order), 'XTickLabel', cellstr(num2str(order)));
    set(gca,'YTick',1:length(order), 'YTickLabel', cellstr(num2str(order)));
    xlabel('Predicted'); ylabel('True');
    title(titleStr);
    for ii=1:size(mat_counts,1)
        for jj=1:size(mat_counts,2)
            text(jj, ii, sprintf('%.2f', mat_counts(ii,jj)), 'HorizontalAlignment','center', 'Color','w', 'FontWeight','bold');
        end
    end
    saveas(fig, fname);
    close(fig);
end

figure(2);
hold on;

passages = [1 3 5 10 15 20 30 50];
colors = lines(length(passages));

for k = 1:length(passages)

    n = passages(k);

    X = clone_size_history{n}; % Have to load clone size history
    X = X(X>0);

    % Relative abundances
    f = X/sum(X);
    fc = mean(f); 
    % PDF
    
    [pdf_vals,edges] = histcounts(f,30,'Normalization','pdf');
    centers = (edges(1:end-1)+edges(2:end))/2;
    % Remove zero PDF bins
    idx = pdf_vals > 0;
    x = centers/fc;
    y = pdf_vals*fc;
    plot(x(idx),log(y(idx)),...
        'o','Color',colors(k,:),...
        'LineWidth',1.5,...
        'MarkerSize',5,...
        'DisplayName',['n = ',num2str(n)]);
end

xlabel('Scaled Relative abundance');
ylabel('Scaled Probability');
legend('Location','best');

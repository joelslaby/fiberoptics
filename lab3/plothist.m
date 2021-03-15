function [] = plothist(data, bins)
    
    mu = mean(data);
    sigma = std(data);

    histogram(data, bins, 'Normalization','pdf'); hold on;

    y = (mu-4*sigma):abs(mu)/1e4:(mu+4*sigma);
    f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
    plot(y,f,'LineWidth',2, 'color', 'r');
    xlim([(mu-4*sigma), (mu+4*sigma)]);

end
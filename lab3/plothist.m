function [] = plothist(data, bins, mean, std)

    histogram(data, bins, 'Normalization','pdf'); hold on;

    y = (mean-4*std):abs(mean)/1e4:(mean+4*std);
    mu = mean;
    sigma = std;
    f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
    plot(y,f,'LineWidth',2, 'color', 'r');
    xlim([(mean-4*std), (mean+4*std)])

end
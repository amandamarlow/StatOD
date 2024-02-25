function plotCovEllipse(cov)
        %# substract mean
        Mu = zeros(2,1);
%         X0 = bsxfun(@minus, X, Mu);

        STD = 3;                     %# 2 standard deviations
        conf = 2*normcdf(STD)-1;     %# covers around 95% of population
        scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions

%         Cov = cov(X0) * scale;
        Cov = cov * scale;
        [V, D] = eig(Cov);
        [D, order] = sort(diag(D), 'descend');
        D = diag(D);
        V = V(:, order);

        t = linspace(0,2*pi,100);
        e = [cos(t) ; sin(t)];        %# unit circle
        VV = V*sqrt(D);               %# scale eigenvectors
        e = bsxfun(@plus, VV*e, Mu); %#' project circle back to orig space

        %# plot cov and major/minor axes
        plot(e(1,:), e(2,:), '--', 'LineWidth', 1.5);
end


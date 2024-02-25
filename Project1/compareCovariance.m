function [] = compareCovariance(P1, P2, axLabels, legLabels, title, pxlim)
%COMPARECOVARIANCE Summary of this function goes here
%   Detailed explanation goes here
% figure
subplot(1,3,1)
sgtitle(title, 'Interpreter', 'latex')
hold on
plotCovEllipse(P1(1:2,1:2,end))
plotCovEllipse(P2(1:2,1:2,end))
grid on
legend(legLabels);
xlabel(axLabels(1))
ylabel(axLabels (2))
xlim([-pxlim pxlim])
axis equal
% hold off
%---------------------------------
subplot(1,3,2)
hold on
plotCovEllipse(P1([1,3],[1,3],end))
plotCovEllipse(P2([1,3],[1,3],end))
grid on
legend(legLabels);
xlabel(axLabels(1))
ylabel(axLabels(3))
xlim([-pxlim pxlim])
axis equal
% hold off
%---------------------------------
subplot(1,3,3)
hold on
plotCovEllipse(P1([2,3],[2,3],end))
plotCovEllipse(P2([2,3],[2,3],end))
grid on
legend(legLabels);
xlabel(axLabels(2))
ylabel(axLabels(3))
xlim([-pxlim pxlim])
axis equal
% hold off
end


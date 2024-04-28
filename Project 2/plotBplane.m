function [] = plotBplane(BdotVec_vec, P_B_vec, legendStrings, plottitle, trueBdotVec, earthRad)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure
hold on
for i = 1:size(BdotVec_vec,2) 
    ellipse = plotBplaneCovEllipse(P_B_vec(:,:,i), BdotVec_vec(:,i));
    center = scatter(BdotVec_vec(1,i),BdotVec_vec(2,i), '*');
    center.SeriesIndex = ellipse.SeriesIndex;
end
if isempty(trueBdotVec)
    if earthRad > 0
        colorE = [0.4660 0.6740 0.1880];
        circle(0,0,earthRad, colorE)
        legend([legendStrings, "Earth"], 'Interpreter', 'latex', 'Location','best')
    else
        legend(legendStrings, 'Interpreter', 'latex', 'Location','best')
    end
else
    scatter(trueBdotVec(1),trueBdotVec(2), '*')
    if earthRad > 0
        colorE = [0.4660 0.6740 0.1880];
        circle(0,0,earthRad, colorE)
        legend([legendStrings, "True Target", "Earth"], 'Interpreter', 'latex', 'Location','best')
    else
        legend([legendStrings, "True Target"], 'Interpreter', 'latex', 'Location','best')
    end
end
axis equal
grid on
xlabel("B$\cdot$T", "Interpreter","latex")
ylabel("B$\cdot$R", "Interpreter","latex")
title(plottitle, "Interpreter","latex")
end


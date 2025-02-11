clc
clear
algorithms = {'1A','2A','3A','4A'};
filename = 'SOS'; run = 31;
algorithmsNumber = length(algorithms); functionsNumber = 1; experimentNumber = 1; 
solutionR = zeros(algorithmsNumber, functionsNumber * experimentNumber, run);
solution = zeros(algorithmsNumber * experimentNumber, functionsNumber, run);

for i = 1 : algorithmsNumber
    solutionR(i,:,:) = xlsread(filename, char(algorithms(i)));
end


 array = [3 7 11 23];
for fn=1:functionsNumber
    values = [];
    
    for m = 1  :algorithmsNumber
       % values = [values ;  squeeze(solutionR(m, fn, :))'];
       values(m, :) = solutionR(m, fn, :);
    end
%    boxplot(values','Notch','off','Labels',{'mu = 3','mu = 6', 'mu = 2', 'mu = 64'}, 'OutlierSize',1);
%     k = fn;
%     if k >= 2
%         k = k +1;
%     end


    boxplot(values', {'1A','2A','3A','4A'});
%         title(strcat('F',int2str(k),{' '},siralar(k)));

    xlabel('Algorithms');
    ylabel('Function Error Value');
    saveas(gcf,strcat('boxPlot/',int2str(fn),'.png'));
    
    
    
    
    
end
    
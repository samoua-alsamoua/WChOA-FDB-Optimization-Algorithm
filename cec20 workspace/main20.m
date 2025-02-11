clear;
clc;
%algorithms = {'FDB_WChimp_CASE1'};
algorithms = {'WChOA','FDB_WChimp_CASE1','FDB_WChimp_CASE6','FDB_WChimp_CASE15','FDB_WChimp_CASE27','FDB_WChimp_CASE28','FDB_WChimp_CASE40','FDB_WChimp_CASE53','FDB_WChimp_CASE57'};
dimension = 10; 
run = 1; 
maxFE = 1*dimension;
filename = 'result'; 
functionsNumber = 10;
solution = zeros(functionsNumber, run);
globalMins = [100, 1100, 700, 1900, 1700, 1600, 2100, 2200, 2400, 2500];
paths;
cec20so = str2func('cec20_func_so'); 
for ii = 1 : length(algorithms)
    disp(algorithms(ii));
    algorithm = str2func(char(algorithms(ii)));
    for i = 1 : functionsNumber
        disp(i);
        for j = 1 : run 
            [bestSolution, bestFitness, iteration] = algorithm(cec20so, dimension, maxFE, i);
            solution(i, j) = bestFitness - globalMins(i);
   
        end
    end
    xlswrite(strcat(filename, '-d=', num2str(dimension), '.xlsx'), solution, func2str(algorithm));
    eD = strcat(func2str(algorithm), '-Bitti :)');
    disp(eD);
end
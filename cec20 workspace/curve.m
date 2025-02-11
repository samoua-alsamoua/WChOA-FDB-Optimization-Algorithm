paths; 
[cec20_func_so] = configuration();
%algorithms = {'sfs_case2', 'lshade_spacma', 'lshade', 'lshade_cnepsin', 'efo', 'mrfo'};
algorithms = {'sos'};
dimension = 30; dim = 1;
% phase = [999, 699, 1399];
% index = [1000, 1700, 3400];
phase = [9999, 6999, 13999];
index = [10000, 17000, 34000];
maxIteration = 10000*dimension;
fitnessA = zeros(maxIteration, length(algorithms));
functionNumber = 1;

for i = 1 : length(algorithms)
    disp(algorithms(i));
    algorithm = str2func(strcat('curve_', char(algorithms(i))));
    fitnessA(:,i) = algorithm(cec20_func_so, dimension, maxIteration, functionNumber) - functionNumber * 100;
end

tA=1:1:maxIteration;
fitness = downsample(fitnessA, index(dim), phase(dim));
fitness = [fitnessA(1, :); fitness];
t = downsample(tA, index(dim), phase(dim));
t = [1, t];

figure('Position',[0 0 640 480])
loglog(t,fitness(:,1),'-k*',':md');
% plot(t,fitness(:,1),'-k*',t,fitness(:,2),':r+',t,fitness(:,3),':bo',t,fitness(:,4),':gx',t,fitness(:,5),':bs',t,fitness(:,6),':gd',t,fitness(:,7),':r.',t,fitness(:,8),'--g+',t,fitness(:,9),'--ro',t,fitness(:,10),'--bx',t,fitness(:,11),'--rd',t,fitness(:,12),'--gs',t,fitness(:,13),'--b.');
% loglog(t,fitness(:,1),'-k*',t,fitness(:,2),':ro',t,fitness(:,3),':bd',t,fitness(:,4),':gs',t,fitness(:,5),':bo',t,fitness(:,6),':gd',t,fitness(:,7),':rs');
% ylim([0 0.1]);
xlabel('Number of function evaluations');
ylabel('Function error values');
legend({'FDBSFS', 'LSHADE\_SPACMA', 'LSHADE', 'LSHADE\_CNEPSIN', 'EFO', 'MRFO'},'Location','southwest');
% legend('SFS', 'CASE1', 'CASE2', 'CASE3');
disp('Bitti');

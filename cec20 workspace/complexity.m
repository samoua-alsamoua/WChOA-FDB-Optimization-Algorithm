paths; 
algorithms = {'sfs_case2', 'lshade_spacma', 'lshade', 'lshade_cnepsin', 'efo', 'mrfo'};
algorithmsNumber = length(algorithms); cec17 = str2func('cec17_func');
dimensions = [30, 50, 100];
[~, dSize] = size(dimensions);
t1 = zeros(1, dSize); 
t2 = zeros(1, 5);

tic;
x = 0.55;
for i = 1 : 1000000
    x = x + x; 
    x = x / 2; 
    x = x * x; 
    x = sqrt(x); 
    x = log(x); 
    x = exp(x); 
    x = x / (x + 2);
end
t0 = toc;
disp(t0);

for ii = 1 : dSize
    dimension = dimensions(ii);
    disp(dimension);
    y = zeros(1, dimension);
    tic;
    for i = 1 : 200000
        fx = feval(cec17, y', 18);
    end
    t1(ii) = toc;
    disp(t1(ii));
    
    for i = 1 : algorithmsNumber
        algorithm = str2func(char(algorithms(i)));
        for times = 1 : 5
            tic;
            algorithm(cec17, dimension, 200000, 18);
            t2(times) = toc;
        end
        t2r = (mean(t2)-t1(ii))/t0;
        disp(algorithms(i));
        disp(t2r);
    end
end

disp('Bitti');
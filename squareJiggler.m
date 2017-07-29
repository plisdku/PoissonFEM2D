% Generate a pretty-random mesh

n_side = 10;

xx = zeros(n_side);
yy = zeros(n_side);

alpha = 0.25;
beta = 1 - 0.5*alpha;

for cc = 1:n_side-1
    
    xx(1:cc, cc+1) = xx(1:cc, cc) + beta + alpha*rand(cc,1);
    yy(1:cc, cc+1) = yy(1:cc, cc) - 0.5*alpha + alpha*rand(cc,1);
    
    %disp(yy)
    %pause
    
    xx(cc+1,1:cc) = xx(cc,1:cc) - 0.5*alpha + alpha*rand(1,cc);
    yy(cc+1,1:cc) = yy(cc,1:cc) + beta + alpha*rand(1,cc);
    
    %disp(yy)
    %pause
    
    xx(cc+1,cc+1) = mean(xx(1:cc,cc+1)) - 0.5*alpha + alpha*rand();
    yy(cc+1,cc+1) = mean(yy(cc+1,1:cc)) - 0.5*alpha + alpha*rand();
    
    %disp(yy)
    %pause
    
end
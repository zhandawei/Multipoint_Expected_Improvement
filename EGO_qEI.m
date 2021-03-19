% EGO with qEI infill criterion
clearvars;clc;close all;
addpath('dace');
% problem settings
fun_name = 'GoldPrice';
num_vari = 2;
lower_bound = [-2,-2]; 
upper_bound =  [2,2]; 
% algorithm settings
num_initial = 10*num_vari;
num_q = 2;
max_iteration = 40;
% number of current iteration
iteration = 1;
% generate random samples
sample_x = lhsdesign(num_initial, num_vari,'criterion','maximin','iterations',1000).*(upper_bound - lower_bound) + lower_bound;
sample_y = feval(fun_name, sample_x);
evaluation =  size(sample_x,1);
% best objectives in each generation
fmin_record = zeros(max_iteration + 1,1);
fmin = min(sample_y);
fmin_record(iteration,:) = fmin;
% print the iteration information
fprintf('%dEI on %d-D %s, iteration: %d, evaluation: %d, best: %0.4g, maximum qEI: %0.8g\n',num_q,num_vari,fun_name,iteration-1,evaluation,fmin,0);
% the evoluation of the population
while iteration <= max_iteration
    % build the Kriging model
    kriging_model = dacefit(sample_x,sample_y,'regpoly0','corrgauss',1*ones(1,num_vari),0.01*ones(1,num_vari),100*ones(1,num_vari));
    % using GA to find the optimum of qEI
    options = gaoptimset('PopulationSize',20,'Generations',50,'StallGenLimit',50,'Display','off','Vectorized','off');
    infill_criterion = @(x)-Infill_qEI(x,kriging_model,fmin);
    [best_x,max_EI]= ga(infill_criterion,num_vari*num_q,[],[],[],[],repmat(lower_bound,1,num_q),repmat(upper_bound,1,num_q),[],[],options);
    infill_x = reshape(best_x,num_vari,[])';
    % evaluate the infill samples in parallel
    infill_y = feval(fun_name, infill_x);
    % update database
    sample_x = [sample_x;infill_x];
    sample_y = [sample_y;infill_y];
    % update the evaluation number of generation number
    iteration = iteration + 1;
    evaluation = evaluation + size(infill_x,1);
    [fmin,index] = min(sample_y);
    fmin_record(iteration,:) = fmin;
    % print the iteration information
    fprintf('%dEI on %d-D %s, iteration: %d, evaluation: %d, best: %0.4g, maximum qEI: %0.8g\n',num_q,num_vari,fun_name,iteration-1,evaluation,fmin,-max_EI)
end



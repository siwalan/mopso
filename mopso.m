%% Multiple Objective Particle Swarm Optimization
% Modified according to known best practice.
% Application of Computational Intelligence - As Taught Professor I. T. Yang
% National Taiwan University of Science and Technology
% Created by : Danny Gho

clc; clear;
w_init = 1;c1 = 2; c2 =2;
pop_size = 100; iteration = 10;
% Boundary Limit for the Particles, Please change it to match your condition.
ub = [5 3];
lb = [0 0];
vub = [0.5 0.5];
vlb = [-0.5 -0.5];

n_var = length(ub);

mo_fit = 'NNDI';
%% Particle Initialization
x = rand(pop_size, n_var) .* (ub-lb) + lb;

%% Initial Fitness Evaluation
fv = [];
for pop_iter=1:pop_size
 fv(pop_iter,:) = objFunction(x(pop_iter,:));
end

pbestfv = fv;
pbest = x;

%%% Check for Domination
[pdominant, pdominantfv] = sortDomination(x, fv, pop_size)

% Determine the Fitness of the dominant
% Nearest Neighbor Density Estimator
pdominantfv_NNDE = distNNDE(pdominantfv);

[value, index] = max(pdominantfv_NNDE);

gbestNNDE = value;
gbestfv = pdominantfv(index,:);
gbest = pdominant(index,:);

v = (rand(pop_size, n_var)) .* (vub-vlb) + vlb;
%% Population Iteration
for iterator = 1:iteration
w = w_init - (0.5/iteration * iterator);

    for pop_iter=1:pop_size

        %Determine the Velocity
        v(pop_iter,:)  = w * v(pop_iter,:)  + rand(1,n_var) .* c1 .* (pbest(pop_iter,:)-x(pop_iter,:)) + rand(1,n_var)  .* c2 .* (gbest - x(pop_iter,:));

        % Calculate the new Position
        x(pop_iter,:) = x(pop_iter,:) + v(pop_iter,:) ;

        % Fixing the Boundary
        bindex_up = x(pop_iter,:) > ub;
        bindex_down = x(pop_iter,:) < lb;

        x(pop_iter,bindex_up)=ub(bindex_up);
        x(pop_iter,bindex_down)=lb(bindex_down);

        % Calculate the New Fitness Value
        fv(pop_iter,:) = objFunction(x(pop_iter,:)); 

        % Damping the Velocity
        vbindex_up = v(pop_iter,:) > vub;
        vbindex_down = v(pop_iter,:) < vlb;

        v(pop_iter,vbindex_up)=vub(vbindex_up);
        v(pop_iter,vbindex_down)=vlb(vbindex_down);

        dominated = checkDomination(fv(pop_iter,:), pbestfv(pop_iter));
        valid = checkValid(x(pop_iter,:),lb,ub);
        if (dominated == 0 & valid == 1)
             pbestfv(pop_iter,:) = fv(pop_iter,:);
             pbest(pop_iter,:) = x(pop_iter,:);
        end
        
        
    end
    
    % Expand the Domination Matrix
    [pdominant, pdominantfv] = sortDomination([pdominant; x], [pdominantfv; fv], pop_size+size(pdominantfv,1))

    %% Determine the GBest
    pdominantfv_NNDE = distNNDE(pdominantfv);
    [value, index] = max(pdominantfv_NNDE);

    if value > gbestNNDE
        gbestfv = pdominantfv(index,:);
        gbest = pdominant(index,:);
    end
end

if (nvar == 2)
   scatter(pdominantfv(:,1), pdominantfv(:,2)) 
end

function fitness = objFunction(var)
    %%% Test Case from http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.46.8661
    fitness(1) = 4*var(1)^2 + 4*var(2)^2;
    fitness(2) = (var(1) -5)^2 + (var(2)-5)^2;
end


function valid = checkValid(particle, lb, ub)
    % Check for if the particle range is under UB and over LB
    break_boundary = sum(particle>=ub) + sum(particle<=lb);
    
    % TO DO
    % Handle Constraint!!! 
    
    if (break_boundary == 0)
        valid = 1;
    else
        valid = 0;
    end
end

function [dominated, dominating] = checkDomination(particle,particles)
    dominated = sum(sum((particles<=particle),2) == length(particle));
    dominating = sum(sum((particles>=particle),2) == length(particle));
end

function [pdominant, pdominantfv] = sortDomination(x, fv, pop_size)
    pdominant = [];
    pdominantfv = [];
    for pop_iter=1:pop_size
        fv_iter = fv;
        fv_iter(pop_iter,:) = [];
        dominated = checkDomination(fv(pop_iter,:), fv_iter);
        if (dominated == 0)
           pdominantfv = [pdominantfv; fv(pop_iter,:)];
           pdominant = [pdominant; x(pop_iter,:)];
        end
    end
end

function pdominantfv_NNDE = distNNDE(pdominantfv)
    pdominantfv_NNDE = []
    
    for iterator=1:size(pdominantfv,1)
        pdominantfv_NNDE = [pdominantfv_NNDE; sum(mink(sqrt(sum((pdominantfv(iterator,:)-pdominantfv(:,:)).^2,2)),3))];    
    end
end

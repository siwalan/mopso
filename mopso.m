%% Multiple Objective Particle Swarm Optimization
% Modified according to known best practice.
% Application of Computational Intelligence - As Taught Professor I. T. Yang
% National Taiwan University of Science and Technology
% Created by : Danny Gho

clc; clear;
w_init = 1;c1 = 2; c2 =2;
pop_size = 1000; iteration = 200;
ub = [5 3];
lb = [0 0];
vub = [0.5 0.5];
vlb = [-0.5 -0.5];

n_var = length(ub);

figure_iterator = 1;
%% Particle Initialization
x = rand(pop_size, n_var) .* (ub-lb) + lb;

%% Initial Fitness Evaluation
fv = [];
for pop_iter=1:pop_size
 fv(pop_iter,:) = objFunction(x(pop_iter,:));
end

%%% We assume that all the particle is... dominating agaisnt it previous
%%% self.
pbestfv = fv;
pbest = x;
fv_history(:,:,1) = fv;
x_history(:,:,1) = x;

%%% Check for Domination
pdominant = [];
pdominantfv = [];
for pop_iter=1:pop_size
    pbestfv_trial = pbestfv;
    pbestfv_trial(pop_iter,:) = [];
    dominated = checkDomination(pbestfv(pop_iter,:), pbestfv_trial);
    if (dominated == 0)
       pdominantfv = [pdominantfv; pbestfv(pop_iter,:)];
       pdominant = [pdominant; x(pop_iter,:)];
    end
end

pdominantfv_history = [];
pdominant_history = [];

%% Determine the Fitness of the dominant
%%% Nearest Neighbor Density Estimator
pdominantfv_NNDE = [];
for iterator=1:size(pdominantfv,1)
    pdominantfv_NNDE = [pdominantfv_NNDE; sum(mink(sqrt(sum((pdominantfv(iterator,:)-pdominantfv(:,:)).^2,2)),3))];    
end
[value, index] = max(pdominantfv_NNDE);

gbestNNDE = value;
gbestfv = pdominantfv(index,:);
gbest = pdominant(index,:);

pdominantfv_history = pdominantfv;
pdominant_history = pdominant;

v = (rand(pop_size, n_var)) .* (vub-vlb) + vlb;
%% Population Iteration
for iterator = 1:iteration
w = w_init - (0.5/iteration * iteration);

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
        if (dominated == 0 & valid == 1)
             pbestfv(pop_iter,:) = fv(pop_iter,:);
             pbest(pop_iter,:) = x(pop_iter,:);
        end
        
        fv_history(:,:,1+iterator) = fv;
        x_history(:,:,1+iterator) = x;
        
    end
    
    %% Expand the Domination Matrix
    for dominant_loop_iterator=1:pop_size
        pbestfv_trial = pbestfv;
        pbestfv_trial(dominant_loop_iterator,:) = [];
        dominated = checkDomination(pbestfv(dominant_loop_iterator,:), pbestfv_trial);
        if (dominated == 0)
           pdominantfv = [pdominantfv; pbestfv(dominant_loop_iterator,:)];
           pdominant = [pdominant; pbest(dominant_loop_iterator,:)];
        end
    end
    
    %% Remove Duplicate Dominanating Resut
    [B,I] =unique(pdominant,'rows');
    pdominant = B;
    pdominantfv = pdominantfv(I,:);
    
    %% Kill Non Dominating Result
    dominated_rank = [];
        for dominant_loop_iterator=1:size(pdominantfv,1)
            pdominantfv_trial = pdominantfv;
            pdominantfv_trial(dominant_loop_iterator,:) = [];
            dominated = checkDomination(pdominantfv(dominant_loop_iterator,:), pdominantfv_trial);
            dominated_rank = [dominated_rank; dominated];
        end
        
    I = find(dominated_rank == 0);
    pdominantfv = pdominantfv(I,:);
    pdominant = pdominant(I,:);

    %% Determine the GBest
    pdominantfv_NNDE =[];
    for NNDE_iterator=1:size(pdominantfv,1)
        pdominantfv_NNDE = [pdominantfv_NNDE; sum(mink(sqrt(sum((pdominantfv(NNDE_iterator,:)-pdominantfv(:,:)).^2,2)),3))];    
    end

    [value, index] = max(pdominantfv_NNDE);

    if value > gbestNNDE
        gbestfv = pdominantfv(index,:);
        gbest = pdominant(index,:);
    end
    
    pdominantfv_history = [pdominantfv_history; 0 0;pdominantfv];
    pdominant_history = [pdominant_history; 0 0;pdominant];
end


function fitness = objFunction(var)
    %%% Test Case from http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.46.8661
    fitness(1) = 4*var(1)^2 + 4*var(2)^2;
    fitness(2) = (var(1) -5)^2 + (var(2)-5)^2;
end


function valid = checkValid(particle, lb, ub)
   break_boundary = sum(particle>=ub) + sum(particle<=lb);
   break_boundary = sum(particle>=ub) + sum(particle<=lb);
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


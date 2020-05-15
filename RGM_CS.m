%_________________________________________________________________________%
%  RGM-MOA framework with cuckoo search (RGM-CS) source codes             %
%                                                                         %
%  Developed in MATLAB                                      %
%                                                                         %
%  Author and programmer: Zhesen Cui                              %
%                                                                         %
%         E-Mail: cuizhesen@gmail.com
%_________________________________________________________________________%
function main()
clc;clear all; close all;

E0 = 0.001;
Itermax = 100;
narvs = 1;
Nests = 50;     % number of nest
lambda = 1.5;   % power law index
pa=0.25;        % Discovery rate of alien eggs/solutions
alpha = 1;      % Step Length
Tolerance = 0.9;   % Minimum Fitness

dim=1;
lb=0;
ub=1;
load raw_data
[M,N] = size(raw_data);
result_rgm=zeros(1,2);
PRE_all=zeros(6,1);
for sample=1:M
    input=raw_data(sample,1:5);
    output=raw_data(sample,6);
    x=initialization(Nests,1,1,0);
    weights = x';
    for i=1:Nests
        [PRE MRE weights(i,:)]=RollingGM(weights(i,:),input,output);
        fitness1(i)=MRE;
    end
    cuckoo = lb+(ub-lb).*rand(1,dim);   % Cuckoo initialized
    minfit = 0;
    iter = 0;
    while (iter < Itermax) && (minfit == 0)
        iter = iter +1;
        cuckoo = cuckoo + alpha.*levy_cs(1,dim,lambda); % Levy flight
        [PRE MRE]=RollingGM(cuckoo,input,output);
        cuckoo_fitness = MRE;
        j = randi([1 Nests]);  % Selecting random nest
        if cuckoo_fitness < fitness1(j)
            weights(j,:) = cuckoo;
            fitness1(j) = cuckoo_fitness;
        end
        % Arranging in descending order
        temp = [(fitness1)' weights];
        temp_sort = sortrows(temp,1);
        new_fitness = (temp_sort(:,1))';
        new_weights = temp_sort(:,2:end);
        % Removing worst nests and building new ones
        remove = int16(pa*Nests);
        weights = new_weights(1:(Nests-double(remove)),:);
        join_weights = lb+(ub-lb).*rand(remove,dim);
        weights = [weights;join_weights];    % updated weights
        
        if new_fitness(1) > Tolerance
            minfit = 1;
        end
        
        for i=1:Nests
            [PRE MRE weights(i,:)]=RollingGM(weights(i,:),input,output);
            fitness1(i)=MRE;
        end
        cuckoo = weights(end,:);
        Best_score = new_fitness(:,1);
    end
    best_soln = weights(1,:);
    result_rgm=[result_rgm;best_soln new_fitness(:,1)];
end

result_rgm=result_rgm(2:end,1);
for i=1:M
    input=raw_data(i,1:5);
    output=raw_data(i,6);
    lambda=result_rgm(i);
    [PRE MRE lambda]=RollingGM(lambda,input,output);
    PRE_all=[PRE_all,PRE'];
end

forecasting_results=PRE_all(6,2:end);
fprintf('forecasting values for tourism demands from 2013 to 2018\n')
forecasting_results
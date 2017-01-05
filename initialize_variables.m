function chromosome = initialize_variables(chromosome,varargin)
global CostFunction pop nVar VarMin VarMax numOfObj TestProblem dynamic;
min = VarMin;
max = VarMax;

% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover and mutation only the decision variables are used
% while for selection, only the objective variable are utilized.

K = numOfObj + nVar;
%% Initialize each chromosome
% For each chromosome perform the following (N is the population size)
index = randperm(pop,round(varargin{1}*pop));
for i = index
    for j = 1 : nVar
        chromosome(i,j) = min(j) + (max(j) - min(j))*rand(1);
    end
end
for i = 1:pop
    chromosome(i,nVar + 1: K) = CostFunction(chromosome(i,1:nVar));
end
clear index i;
end

function DNSGA_II_B()
global pop colors itrCounter TestProblem step window dynamic CostFunction nVar VarMin VarMax numOfObj;
% Original algorithm NSGA-II was developed by researchers in Kanpur Genetic
% Algorithm Labarotary and kindly visit their website for more information
% http://www.iitk.ac.in/kangal/
colors = {'bo','go','ro','co','mo','ko','bv','gv','rv','cv','mv','kv','bs','gs','rs','cs','ms','ks'};
step = 10;
window = 20;
itrCounter = 1;
TestProblem = 37;
pop = 100;
pool = round(pop/2);
tour = 2;
mu = 10;
mum = 20;
mumrate = 0.2;
maxIt = round(step*window);
[numOfObj, nVar, VarMin, VarMax] = objective_description_function();
chromosome = [];
chromosome = initialize_variables(chromosome, 1);
chromosome = non_domination_sort_mod(chromosome, numOfObj, nVar);

%%%%%%%%%%
for itrCounter = 1 : maxIt
    parent_chromosome = tournament_selection(chromosome, pool, tour);
    offspring_chromosome = ...
        genetic_operator(parent_chromosome, ...
        numOfObj, nVar, mu, mum, VarMin, VarMax, itrCounter);
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    clear temp
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : numOfObj+nVar) = ...
        offspring_chromosome;
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, numOfObj, nVar);
    chromosome = replace_chromosome(intermediate_chromosome, numOfObj, nVar, pop);
    if ~mod(itrCounter,10)
        disp([num2str(itrCounter),' generations completed, chromesome number ',num2str(size(chromosome,1))])
    end
    population2pic(chromosome);
    if mod(itrCounter,window) == 0 && dynamic == 1
        chromosome = mutate(chromosome,mumrate,mum);
        chromosome = non_domination_sort_mod(chromosome, numOfObj, nVar);
    end
end
population2pic(chromosome,chromosome);
%% Visualize

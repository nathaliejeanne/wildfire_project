% Initialisation ----------------------------------------------------------
clc; clear;

setenv("GUROBI_VERSION", "110")
import casadi.*

N = 6; % height/width of the forest
T = 8; % number of the timesteps, i.e. t=1 : the fire starts, t=2 : next timestep, etc.
max_number_of_firefighters = 16; % max number of firefighters
num_poss_starting_fire = 5; % number of possibilities for starting fires

% Different possibilities for a starting fire
i_0 = [1 2 3 5 4];
j_0 = [2 4 3 2 2];

% Different possibilities for a second starting fire
i_0_2 = [2 6 3 5 4];
j_0_2 = [2 6 4 5 4];


% Solutions of the IPs ----------------------------------------------------

number_of_unburned_area = zeros(num_poss_starting_fire,max_number_of_firefighters);
for l = 1:num_poss_starting_fire % go over all possible starting fires
    for K = 1:max_number_of_firefighters % go over all possible max sizes of firefighters
        [x_1, x_2] = fire_problem(N,T,K,i_0(l),j_0(l),i_0_2(l),j_0_2(l),0.01); % solve the IP
        number_of_unburned_area(l,K) = N*N - nnz(x_1{T}); % number of unburned area if K firefighters can be used
    end
end


% Plot --------------------------------------------------------------------
for l = 1:num_poss_starting_fire
    plot(1:K, number_of_unburned_area(l,:)); hold on;
end
hold off;

set(gca, 'FontSize', 16);  % Axis ticks and labels
set(findall(gca, 'Type', 'text'), 'FontSize', 16);  % All text objects
xlim([1 max_number_of_firefighters])
legend({'(1,2) \& (2,2)', '(2,4) \& (6,6)', '(3,3) \& (3,4)', ...
       '(5,2) \& (5,5)', '(4,2) \& (4,4)'}, ...
       'Interpreter', 'latex', 'FontSize', 16, 'Location', 'southeast');
xlabel('$K$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('number of unburned areas', 'Interpreter', 'latex', 'FontSize', 16);


% figure_name = ['two_starting_fires_3.pdf'];
% exportgraphics(gcf, figure_name, 'ContentType', 'vector');
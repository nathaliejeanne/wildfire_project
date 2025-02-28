% INITIALISATION ----------------------------------------------------------
clc; clear;

setenv("GUROBI_VERSION", "110")
import casadi.*

N = 8; % height/width of the forest
T = 8; % number of the timesteps, i.e. t=1 : the fire starts, t=2 : next timestep, etc.

% The fire starts in this area
i_0 = 1;
j_0 = 2;

% A second fire starts in this area (in case of just one fire, write -1)
i_0_2 = 8;
j_0_2 = 7;

% number of firefighters
K = 4;


% SOLUTION OF THE IP ------------------------------------------------------
[grid_fire, grid_firefighters] = fire_problem(N,T,K,i_0,j_0,i_0_2,j_0_2,0.01);


% number of burned and unburned area
number_burned_area = nnz(grid_fire{T})
number_unburned_area = N*N - nnz(grid_fire{T})



% PLOT --------------------------------------------------------------------
number = zeros(N,N);
figure;
for t = 1:T

    grid = grid_fire{t};

    % green = 0 (no fire), red = 1 (fire)
    imagesc(grid);
    colormap([0 1 0; 1 0 0]);
    colorbar off;

    axis equal;
    axis off;
    axis tight;
    % title(['fire spread at the time step t = ', num2str(t), ' and with K = ',num2str(K)],'FontSize', 14);

    % add grid lines
    hold on;
    for i = 0.5:1:N+0.5
        plot([0.5, N+0.5], [i, i], 'k-', 'LineWidth', 0.5);
        plot([i, i], [0.5, N+0.5], 'k-', 'LineWidth', 0.5);
    end

    for i = 1:N
        for j = 1:N

            % Plot of the firefighters
            if grid_firefighters{t}(j,i) == 1
                if number(j,i) == 0
                    text(i, j, num2str(t), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 16, 'FontWeight', 'bold');
                    number(j,i) = t;
                else
                    text(i, j, num2str(number(j,i)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 16, 'FontWeight', 'bold');
                end
            end

            % Numbers of the fire
            if grid_fire{t}(j,i) == 1
                if number(j,i) == 0
                    text(i, j, num2str(t), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 16, 'FontWeight', 'bold');
                    number(j,i) = t;
                else
                    text(i, j, num2str(number(j,i)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 16, 'FontWeight', 'bold');
                end
            end
        end
    end
    
    hold off;
    
    pause(1.0);
end


% figure_name = ['fire_demo_K_4.pdf'];
% exportgraphics(gcf, figure_name, 'ContentType', 'vector');
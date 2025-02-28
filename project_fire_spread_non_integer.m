clc; clear;

import casadi.*

% non-integer optimization variables

N = 6; % height/width of the forest
T = 4; % number of the timesteps, i.e. t=1 : the fire starts, t=2 : next timestep, etc.

% The fire starts in this area
i_0 = 4;
j_0 = 3;

% % A second fire starts in this area
% i_0_2 = 6;
% j_0_2 = 6;

% number of firefighters
K = 3;


opti = casadi.Opti();

x = opti.variable(N*N*T); % if x = 0, the tree is unburned, if x > 0, the tree is burned
u = opti.variable(N*N*T*K); % if u = 0, the tree is unprotected, if u > 0, the tree is protected

% Helper functions to index 3D and 4D variables
index_x = @(i, j, t) (t-1) * (N * N) + (i-1) * N + j;
index_u = @(i, j, t, k) (k-1) * (N * N * T) + (t-1) * (N * N) + (i-1) * N + j;
k_poss = 1:K;

% objective function: minimize the number of burned trees at t = T
opti.minimize( nnz(sum(x(N*N*(T-1)+1:N*N*T))) + 10^(-4)*nnz(sum(u(:))) );

for i = 1:N
    for j = 1:N
        
        % inequality constraint 1 (for time step 1)
        opti.subject_to( {x(index_x(i,j,1)) >= 0, x(index_x(i,j,1)) <= 1} );
        opti.subject_to( {u(index_u(i,j,1,k_poss)) >= 0, u(index_u(i,j,1,k_poss)) <= 1} );

        for t = 2:T
            
            % inequality constraint 2 (now: equality constraint)
            opti.subject_to( x(index_x(i,j,t)) * u(index_u(i,j,t,k_poss)) == 0 );

            % inequality constraint 3 (now: equality constraint)
            opti.subject_to( x(index_x(i,j,t-1)) * ( x(index_x(i,j,t-1))-x(index_x(i,j,t)) ) == 0 );
        
            % inequality constraint 4
            for k = 1:K
               opti.subject_to( u(index_u(i,j,t-1,k)) * (  u(index_u(i,j,t-1,k))-u(index_u(i,j,t,k)) ) >= 0 );
            end
            
            % inequality constraint 6
            neighbours_of_ij = neighbours_of(i,j,N);
            for k = 1:size(neighbours_of_ij,1)
                opti.subject_to( x(index_x(i,j,t)) + sum(u(index_u(i,j,t,k_poss))) - x(index_x(neighbours_of_ij(k,1),neighbours_of_ij(k,2),t-1)) >= 0 );
            end

            % % inequality constraint 7
            % if t > 2
            %     sum_1 = u(index_u(i,j,t,k_poss)) - u(index_u(i,j,t-1,k_poss));
            %     sum_2 = sum(u(index_u(neighbours_of_ij(:,1), neighbours_of_ij(:,2), t-1, k_poss))) - ...
            %     sum(u(index_u(neighbours_of_ij(:,1), neighbours_of_ij(:,2), t-2, k_poss)));
            %     opti.subject_to(sum_1 <= sum_2');
            % end
            
            % inequality constraint 1 (for the time steps t>1)
            opti.subject_to( {x(index_x(i,j,t)) >= 0, x(index_x(i,j,t)) <= 1} );
            opti.subject_to( {u(index_u(i,j,t,k_poss)) >= 0, u(index_u(i,j,t,k_poss)) <= 1} );
        end

        % equality constraint 1
        if (i == i_0) && (j == j_0)
            opti.subject_to( x(index_x(i,j,1)) == 1 );
        % elseif (i == i_0_2) && (j == j_0_2)
        %     opti.subject_to( x(index_x(i,j,1)) == 1 );
        else
            opti.subject_to( x(index_x(i,j,1)) == 0 );
        end

        % equality constraint 2
        opti.subject_to( u(index_u(i,j,1,k_poss)) == 0 );

    end
end

% inequality constraint 5
epsilon = 0.1;
alpha = 50;
for t = 2:T
    for k = 1:K
        summe = 0;
        for i = 1:N
            for j = 1:N
                diff = 0;
                diff = u(index_u(i, j, t, k)) - u(index_u(i, j, t-1, k));
                % summe = summe + diff/(diff+epsilon);
                % summe = summe + (diff^2/(diff^2+epsilon));
                % summe = summe + diff/sqrt(diff^2+epsilon);
                % summe = summe + (1 - exp(-alpha*diff))/(1+exp(-alpha*diff));
                % summe = summe + 1 - exp(-alpha*diff);
                % summe = summe + 1 - exp(-alpha*diff^2);
                % summe = summe + tanh(alpha*diff);
                % summe = summe + 2/pi * atan(alpha*diff);                                
                summe = summe + erf(alpha*diff);
            end
        end
        opti.subject_to( summe <= 1 );
    end
end


opti.solver('ipopt');

tic;
sol = opti.solve();
time = toc


% position of the fire
grid_fire = cell(1, T);
for t = 1:T
    fire_t = sol.value(x((t-1)*N*N+1 : t*N*N));
    grid_fire{t} = reshape(fire_t, N, N)';
end

% position of the firefighters
grid_firefighters = cell(1, T);
for t = 1:T
    firefighters_t = 0;
    for k = 1:K
        firefighters_t = firefighters_t + sol.value(u((k-1)*N*N*T + (t-1)*N*N+1 : (k-1)*N*N*T + t*N*N));
    end
    grid_firefighters{t} = reshape(firefighters_t, N, N)';
end
number = zeros(N,N);

% tolerance
tol = 10^(-1);

% number of burned and unburned area
number_burned_area = nnz(abs(grid_fire{T}) > tol)
number_unburned_area = N*N - number_burned_area


% Plot
figure;
for t = 1:T

    grid = grid_fire{t};

    % green < tol (no fire), red >= tol (fire)
    grid_color = grid > tol;
    imagesc(grid_color);

    colormap([0 1 0; 1 0 0]);
    colorbar off;

    axis equal;
    axis off;
    axis tight;
    title(['fire spread at the time step t = ', num2str(t), ' and with K = ',num2str(K)]);

    % add grid lines
    hold on;
    for i = 0.5:1:N+0.5
        plot([0.5, N+0.5], [i, i], 'k-', 'LineWidth', 0.5);
        plot([i, i], [0.5, N+0.5], 'k-', 'LineWidth', 0.5);
    end

    for i = 1:N
        for j = 1:N

            % Plot of the firefighters
            if grid_firefighters{t}(j,i) >= tol
                if number(j,i) == 0
                    text(i, j, num2str(t), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 12, 'FontWeight', 'bold');
                    number(j,i) = t;
                else
                    text(i, j, num2str(number(j,i)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 12, 'FontWeight', 'bold');
                end
            end

            % Numbers of the fire
            if grid_fire{t}(j,i) >= tol
                if number(j,i) == 0
                    text(i, j, num2str(t), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 12, 'FontWeight', 'bold');
                    number(j,i) = t;
                else
                    text(i, j, num2str(number(j,i)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 12, 'FontWeight', 'bold');
                end
            end
        end
    end
    
    hold off;
    pause(1.0);
end


% Helper function to determine the neighbours of a point (i,j)
function neighbours = neighbours_of(i,j,N)
    if (i==1)&&(j==1)
        neighbours = [1 2; 2 1; 2 2];
    elseif (i==N)&&(j==1)
        neighbours = [N 2 ; N-1 1; N-1 2];
    elseif (i==1)&&(j==N)
        neighbours = [1 N-1; 2 N; 2 N-1];
    elseif (i==N)&&(j==N)
        neighbours = [N N-1; N-1 N; N-1 N-1];
    elseif (i==1)&&(1<j)&&(j<N)
        neighbours = [1 j-1; 1 j+1; 2 j; 2 j-1; 2 j+1];
    elseif (1<i)&&(i<N)&&(j==1)
        neighbours = [i-1 1; i+1 1; i 2; i-1 2; i+1 2];
    elseif (1<i)&&(i<N)&&(j==N)
        neighbours = [i-1 N; i+1 N; i N-1; i-1 N-1; i+1 N-1];
    elseif (i==N)&&(1<j)&&(j<N)
        neighbours = [N j-1; N j+1; N-1 j; N-1 j-1; N-1 j+1];
    else
        neighbours = [i+1 j; i-1 j; i j+1; i j-1; i-1 j-1; i-1 j+1; i+1 j-1; i+1 j+1];
    end
end
% INPUT :

% N            : size of the grid
% T            : number of timesteps
% K            : number of firefighters
% i_0, j_0     : position of the starting fire
% i_0_2, j_0_2 : position of the second starting fire
%                i_0_2 = -1 & j_0_2 = -1 in case there is no second fire
% MIPGap_value : for faster solver time, default value is 1e-4

% OUTPUT :

% x_1_sol      : solution of the IP (fire spread)
% x_2_sol      : solution of the IP (firefighters positions)


function [x_1_sol, x_2_sol] = fire_problem(N,T,K,i_0,j_0,i_0_2,j_0_2,MIPGap_value)

    opti = casadi.Opti('conic');
    
    x_1 = opti.variable(N*N*T); % if x = 1, the tree is burned, if x = 0, the tree is alive (yet)
    x_2 = opti.variable(N*N*T*K); % if u = 1, the tree is or has been defended, if u = 0, not defended (yet)
    opti.set_domain(x_1, 'integer')
    opti.set_domain(x_2, 'integer')
    
    % Helper functions to index 3D and 4D variables
    index_x_1 = @(i, j, t) (t-1) * (N * N) + (i-1) * N + j;
    index_x_2 = @(i, j, t, k) (k-1) * (N * N * T) + (t-1) * (N * N) + (i-1) * N + j;
    k_poss = 1:K;
    
    % objective function: minimize the number of burned trees at t=T
    opti.minimize( sum(x_1(N*N*(T-1)+1:N*N*T)) + 10^(-4)*sum(x_2(:)) );
    
    for i = 1:N
        for j = 1:N
            
            % inequality constraint 1 (for time step 1)
            opti.subject_to( {x_1(index_x_1(i,j,1)) >= 0, x_1(index_x_1(i,j,1)) <= 1} );
            opti.subject_to( {x_2(index_x_2(i,j,1,k_poss)) >= 0, x_2(index_x_2(i,j,1,k_poss)) <= 1} );
    
            for t = 2:T
                
                % inequality constraint 2
                opti.subject_to( x_1(index_x_1(i,j,t)) + x_2(index_x_2(i,j,t,k_poss)) <= 1 );
    
                % inequality constraint 3
                opti.subject_to( x_1(index_x_1(i,j,t)) - x_1(index_x_1(i,j,t-1)) >= 0 );
    
                % inequality constraint 4
                opti.subject_to( x_2(index_x_2(i,j,t,k_poss)) - x_2(index_x_2(i,j,t-1,k_poss)) >= 0 );
    
                % inequality constraint 6
                neighbours_of_ij = neighbours_of(i,j,N);
                for k = 1:size(neighbours_of_ij,1)
                    opti.subject_to( x_1(index_x_1(i,j,t)) + sum(x_2(index_x_2(i,j,t,k_poss))) - x_1(index_x_1(neighbours_of_ij(k,1),neighbours_of_ij(k,2),t-1)) >= 0 );
                end
    
                % inequality constraint 7
                if t > 2
                    sum_1 = x_2(index_x_2(i,j,t,k_poss)) - x_2(index_x_2(i,j,t-1,k_poss));
                    sum_2 = sum(x_2(index_x_2(neighbours_of_ij(:,1), neighbours_of_ij(:,2), t-1, k_poss))) - ...
                    sum(x_2(index_x_2(neighbours_of_ij(:,1), neighbours_of_ij(:,2), t-2, k_poss)));
                    opti.subject_to(sum_1 <= sum_2');
                end
                
                % inequality constraint 1 (for the time steps t>1)
                opti.subject_to( {x_1(index_x_1(i,j,t)) >= 0, x_1(index_x_1(i,j,t)) <= 1} );
                opti.subject_to( {x_2(index_x_2(i,j,t,k_poss)) >= 0, x_2(index_x_2(i,j,t,k_poss)) <= 1} );
            end
    
            % equality constraint 1
            if (i == i_0) && (j == j_0)
                opti.subject_to( x_1(index_x_1(i,j,1)) == 1 );
            elseif (i_0_2 ~= -1) && (j_0_2 ~= -1) && (i == i_0_2) && (j == j_0_2)
                opti.subject_to( x_1(index_x_1(i,j,1)) == 1 );
            else
                opti.subject_to( x_1(index_x_1(i,j,1)) == 0 );
            end
    
            % equality constraint 2
            opti.subject_to( x_2(index_x_2(i,j,1,k_poss)) == 0 );
    
        end
    end
    
    % inequality constraint 5
    for t = 2:T
        for k = 1:K
            opti.subject_to( sum(x_2((k-1)*(N*N*T)+(t-1)*(N*N)+1 : (k-1)*(N*N*T)+(t)*(N*N))) - sum(x_2((k-1)*(N*N*T)+(t-2)*(N*N)+1 : (k-1)*(N*N*T)+(t-1)*(N*N))) <= 1 );
        end
    end

    
    % solve the IP
    opts.gurobi.MIPGap = MIPGap_value;
    opti.solver('gurobi', opts);
    sol = opti.solve();
    
    
    % position of the fire
    grid_fire = cell(1, T);
    for t = 1:T
        fire_t = sol.value(x_1((t-1)*N*N+1 : t*N*N));
        grid_fire{t} = reshape(fire_t, N, N)';
    end
    
    % position of the firefighters
    grid_firefighters = cell(1, T);
    for t = 1:T
        firefighters_t = 0;
        for k = 1:K
            firefighters_t = firefighters_t + sol.value(x_2((k-1)*N*N*T + (t-1)*N*N+1 : (k-1)*N*N*T + t*N*N));
        end
        grid_firefighters{t} = reshape(firefighters_t, N, N)';
    end
    
    % output
    x_1_sol = grid_fire;
    x_2_sol = grid_firefighters;


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

end
% Newton solver for G(d)=F
% required params:
    % argument G should be a function handle for G(d) which returns NxN matrix
    % argument K_R should be a function handle that returns K(d)=Jacibian(G) and R=F-G
    % argument F should be a function handle for F(d). should be of dimension Nx1
% optional params:
    % d_0 is initial guess, defaults to 0 vector. size Nx1
    % max_iterations is max number of iterations. defaults to 300
    % eps_r is required error convergence check. defaults to 1e-5
    % eps_d is required solution convergence check. defaults to 1e-5
% returns vector d
function [solution, I_RC, store_K, store_R, store_D] = newtonSolverF(K_R, d_0)
    global max_iterations eps_r eps_d first_iteration NewtonMethod ID;
    % start time measurment
    tt = cputime;
        
    %% settings
    STORE_K = 1; % set to 1 if we want to record K's iterations
    STORE_R = 1; % set to 1 if we want to record R's iterations
    STORE_D = 1; % set to 1 if we want to record d's iterations

    %% Assert arguments are valid
    if ~isa(K_R, 'function_handle')
        throw(MException('newtonSolver:invalidArguments', 'K_R should be a function handle'));
    end
    N = length(d_0);
    
    %% Begin iterations
    for Niter = 1:max_iterations
        if Niter == 1
            %% first iteration
            first_iteration = 1;
            [next_K,next_R,norm_f] = K_R(d_0);
            next_d = d_0;
            if STORE_K
                calculated_K(:,:,1) = next_K;
            end
            if STORE_R
                calculated_R(:,1) = next_R;
            end
            if STORE_D
                calculated_d(:,1) = d_0;
            end
            
            %% validate right matrix sizes
            if ~isvector(next_R)
                throw(MException('newtonSolver:invalidArguments', 'R should be a vector (Nx1)')); 
            end
        end        
        
        if NewtonMethod(1) || first_iteration
            inverted_K = inv(next_K);
            first_iteration = 0;
        end
        delta_d = mtimes(inverted_K, next_R); % calculate Delta d for unknown nodes
        
        % add delta_d to current d to get next iteration
        for ii=1:length(next_d)
            if ID(ii) > 0
                next_d(ii) = next_d(ii) + delta_d(ID(ii));
            end
        end
        %next_d = next_d + delta_d; % calculate next iteration for d
        
        
        [next_K,next_R,norm_f] = K_R(next_d);   % calculate next iteration for K
        
        %% store results
        if STORE_D
            calculated_d(:,Niter+1) = next_d;
        end
        if STORE_K
            calculated_K(:,:,Niter+1) = next_K;
        end
        if STORE_R
            calculated_R(:,Niter+1) = next_R;
        end
        
        %% Check convergence
        norm_r = norm(next_R)/norm_f;
        norm_d = norm(delta_d)/norm(next_d);
        I_RC(Niter) = 3 - (2*(norm_r < eps_r) + (norm_d < eps_d));
        % ^ sets I_RC code: 0 = full convergence, 1 = error convergence,
        %                   2 = solution convergence, 3 = no convergence
        
        if I_RC(Niter) == 0
           break; 
        end
    end
    
    %% return results
    solution = next_d;
    if STORE_K
        store_K = calculated_K;
    end
    if STORE_R
        store_R = calculated_R;
    end
    if STORE_D
        store_D = calculated_d;
    end
    disp(['Elapsed time: ' num2str(cputime - tt) ' sec']);
end
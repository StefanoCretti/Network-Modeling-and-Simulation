classdef SimulateSSA < BaseSimulation
    % SimulateSSA: simulate system evolution through SSA algorithm
    %     Derived class from BaseSimulation class. 
    %     It uses non-optimized SSA to simulate chemical system evolution. 

    methods
        
        % RUN SSA SIMULATION //////////////////////////////////////////////
        function run_algorithm(obj)
            
            % Create local copies (repeated object property access is slow)
            mat_reag = obj.mod_s.mat_reag;
            mat_evol = obj.mod_s.mat_evol;
            rates    = obj.mod_s.st_rates;
            num_reac = obj.mod_s.num_reac;
            num_mols = obj.mod_s.num_mols;
            init_pop = obj.mod_s.init_pop;
            t_max    = obj.par_s.t_max;
            num_pts  = obj.par_s.n_discr_pt;
            
            % Define some values to decide whether to store the pt 
            curr_pt    = 1;
            step_width = t_max / num_pts;
            curr_thr   = curr_pt*step_width;
            prev_dist  = Inf; 

            % Create simulation results containers
            times    = nan(num_pts, 1);
            dynam    = nan(num_pts, num_mols);

            % Initialize step counter
            num_steps = 0;

            % Store to minimize lookup
            curr_state  = init_pop(1,:);
            curr_time   = 0;
            
            % Simulate until stop condition
            while curr_time < t_max

                % Recompute ALL propensity rates
                prop = rates;
                for r = 1:num_reac
                    k_reac = mat_reag(r,:);
                    for m = 1:num_mols
                        n = curr_state(m);
                        k = k_reac(m);
                        if n < k
                            prop(r) = 0;
                            break
                        elseif k == 1
                            prop(r) = prop(r)*n;
                        else
                            prop(r) = prop(r)*nchoosek(n,k);
                        end
                    end
                end

                tot_prop = sum(prop);
                
                % Define firing reaction
                prop_sum = 0;
                to_fire  = 0;
                fire_thr = tot_prop*rand(1); 
                
                while prop_sum < fire_thr
                    to_fire = to_fire + 1;
                    prop_sum = prop_sum + prop(to_fire);
                end 
                
                % Define time step
                tau = (1/tot_prop)*log(1/rand(1));
                
                % Increase step counter
                num_steps = num_steps+1;
                
                % Update state and time
                prev_time  = curr_time;
                prev_state = curr_state;
                curr_time  = curr_time  + tau;
                curr_state = curr_state + mat_evol(to_fire,:);
                
                % Decide whether to store point
                curr_dist = abs(curr_thr - curr_time);
                
                while (curr_dist > prev_dist) && (curr_pt <= num_pts)

                    % Store previous point 
                    times(curr_pt)    = prev_time ;
                    dynam(curr_pt, :) = prev_state;

                    % Increment point counter and recompute
                    % NOTE: this way if previous simulated point was the 
                    % closest to the new discretized point it is caught
                    curr_pt        = curr_pt + 1;
                    curr_thr       = curr_pt * step_width;
                    prev_dist      = abs(curr_thr - prev_time);
                    curr_dist      = abs(curr_thr - curr_time);
       
                end

                prev_dist = curr_dist;
                
            end % end of simulation
            
            % If t_max is reached, the last point is repeated to pad the
            % number of preallocated fixed time points
            while curr_pt < num_pts + 1
                times(curr_pt)   = prev_time ;
                dynam(curr_pt,:) = prev_state;
                curr_pt = curr_pt + 1;
            end

            % Store simulation results
            obj.res_s.full_tim = times(:,:);
            obj.res_s.full_mat = dynam(:,:);
             
            obj.res_s.sim_stat.num_tot_steps = num_steps;
            obj.res_s.sim_stat.num_ssa_steps = num_steps; 

        end
        %//////////////////////////////////////////////////////////////////

    end

end

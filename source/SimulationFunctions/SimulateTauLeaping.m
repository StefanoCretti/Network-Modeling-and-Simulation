classdef SimulateTauLeaping < BaseSimulation
    % SimulateTauLeaping: simulate system evolution through tau-leaping
    %     Derived class from BaseSimulation class. 
    %     It uses non-optimized tau-leaping algorithm to simulate chemical 
    %     system evolution.


    methods

        % RUN TAU LEAPING SIMULATION //////////////////////////////////////
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

            crit_thr = obj.par_s.n_crit_thr;
            neg_err  = obj.par_s.epsilon;
            
            % Define some values to decide whether to store the pt 
            curr_pt    = 1;
            step_width = t_max / num_pts;
            curr_thr   = curr_pt*step_width;
            prev_dist  = Inf; 

            % Create simulation results containers
            times    = nan(num_pts, 1);
            dynam    = nan(num_pts, num_mols);

            % Store to minimize lookup
            curr_state  = init_pop(1,:);
            curr_time   = 0;

            % Initiate counter of simulation steps
            tau_steps = 0;
            rej_steps = 0;
            
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

                % Recompute ALL critical reactions
                is_crit = false(1,num_reac);
                for r = 1:num_reac
                    if prop(r) > 0
                        lim_firing = inf;
                        evol_coefs = mat_evol(r,:);
                        for m = 1:num_mols
                            mol_coef = abs(evol_coefs(m));
                            if mol_coef > 0
                                num_firing = floor(curr_state(m)/mol_coef); 
                                if num_firing < lim_firing
                                    lim_firing = num_firing;
                                end
                            end
                        end
                        if lim_firing <= crit_thr
                            is_crit(r) = true;
                        end
                    end
                end

                % Code in between hash-lines could be computed only once, 
                % but it is computed at each iteration since the original
                % code does the same (in order to compare results)

                %##########################################################
                % Compute highest stoich coef for each molecular species
                max_coef_m = max(mat_reag);
                % Indexer of reactant species (non zero in at least 1 reac)
                is_react = max_coef_m > 0;
               
                % Compute max order of reactions involving each reactant
                % and maximum coef of the reactant among those reactions
                max_ord_r = nan(1,num_mols);
                max_stoic = nan(1,num_mols);
                for m = 1:num_mols
                    if is_react(m)
                        % Retrieve reactions involving the molecule
                        subset_mat = mat_reag(mat_reag(:,m)>0,:);
                        % Compute reaction orders of all those reactions
                        reac_ords = sum(subset_mat,2);
                        % Compute max ord of any reaction involv. reactant
                        max_ord = max(reac_ords);
                        % Store that value in max_ord_r
                        max_ord_r(m) = max_ord;
                        % Subset to max order reactions for the reactant
                        max_subset = subset_mat(reac_ords == max_ord,:); 
                        % Compute max stoich coef for m among those reacts
                        max_stoic(m) = max(max_subset(:,m));
                    end
                end
                %##########################################################
             
                % Compute gi values
                gi_vals = max_ord_r;
                for m=1:num_mols
                    if is_react(m)
                        max_coef = max_stoic(m);
                        if max_coef > 1
                            add_term = 0;
                            for j=1:(max_coef-1)
                                add_term = add_term + j/(curr_state(m)-j);
                            end
                            add_term = add_term*max_ord_r(m)/max_coef;
                            gi_vals(m) = gi_vals(m) + add_term;
                        end
                    end
                end

                % Compute mi and si^2 values 
                mi_vals = zeros(1,num_mols);
                si_vals = zeros(1,num_mols);
                for m=1:num_mols
                    if is_react(m)
                        m_evol_coefs = mat_evol(:,m);
                        for r=1:num_reac
                            coef_evol = m_evol_coefs(r);
                            if (~is_crit(r)) && (coef_evol ~= 0)
                                r_prop = prop(r);
                                mi_vals(m) = mi_vals(m)+coef_evol*r_prop;
                                si_vals(m) = si_vals(m)+coef_evol^2*r_prop;
                            end
                        end
                    end
                end

                % Compute tau_NC
                tau_NC = inf; 
                for m=1:num_mols
                    if is_react(m)
                        mi = mi_vals(m);
                        numtor = max(neg_err*curr_state(m)/(gi_vals(m)),1);
                        fterm = numtor   / abs(mi);
                        sterm = numtor^2 /(si_vals(m)-abs(mi^2/tot_prop));
                        m_tau = min(fterm,sterm); 
                        if m_tau < tau_NC
                            tau_NC = m_tau;
                        end
                    end
                end       
                
                % Iterate until one value of tau is accepted
                while true 
                    
                    % Condition to break out of the loop
                    accept_step = true; 

                    % Compute number of firings for non-critical reacts
                    k_vals = zeros(1,num_reac);
                    k_vals(~is_crit) = poissrnd(prop(~is_crit)*tau_NC);

                    % Try updating the system with the num of firings
                    % If any negative, divide tau by 2 and restart loop
                    new_state = curr_state + k_vals * mat_evol;
                    for m=1:num_mols
                        if new_state(m) < 0
                            accept_step = false;
                            rej_steps = rej_steps + 1;
                            tau_NC = tau_NC/2;
                            break
                        end
                    end
                    
                    % Proposed update is accepted, save sys evolution
                    if accept_step
                        
                        % Increate t-step counter
                        tau_steps = tau_steps + 1;

                        % Update state and time
                        prev_time  = curr_time;
                        prev_state = curr_state;
                        curr_time  = curr_time  + tau_NC;
                        curr_state = new_state;
                        
                        % Decide whether to store point
                        curr_dist = abs(curr_thr - curr_time);
                        
                        while (curr_dist>prev_dist) && (curr_pt<=num_pts)
            
                            % Store previous point 
                            times(curr_pt)    = prev_time ;
                            dynam(curr_pt, :) = prev_state;
            
                            % Increment point counter and recompute
                            % NOTE: this way if previous simulated pt was  
                            % closer to the new discretized pt it is caught
                            curr_pt        = curr_pt + 1;
                            curr_thr       = curr_pt * step_width;
                            prev_dist      = abs(curr_thr - prev_time);
                            curr_dist      = abs(curr_thr - curr_time);
               
                        end
            
                        prev_dist = curr_dist;
                        
                        % break and recompute critical partition
                        break 

                    end
               
                end % end of tau acceptance/rejection

            end % end of simulation (t_max reached)

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
            
            obj.res_s.sim_stat.num_tau_steps = tau_steps;
            obj.res_s.sim_stat.num_rej_steps = rej_steps;
            obj.res_s.sim_stat.num_tot_steps = tau_steps + rej_steps;

        end
        %//////////////////////////////////////////////////////////////////
    
    end

end

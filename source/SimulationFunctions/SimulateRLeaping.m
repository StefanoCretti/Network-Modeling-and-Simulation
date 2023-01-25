classdef SimulateRLeaping < BaseSimulation
    % SimulateRLeaping: simulate system evolution through R-leaping
    %     Derived class from BaseSimulation class. 
    %     It uses non-optimized R-leaping algorithm to simulate chemical 
    %     system evolution. 
    
    methods
       
        % RUN R LEAPING SIMULATION ////////////////////////////////////////
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

            neg_err  = obj.par_s.epsilon;
            sort_fr  = obj.par_s.reord_freq;
            theta    = obj.par_s.theta;
            
            % Define some values to decide whether to store the pt 
            curr_pt    = 1;
            step_width = t_max / num_pts;
            curr_thr   = curr_pt*step_width;
            prev_dist  = Inf; 

            % Create simulation results containers
            times    = nan(num_pts, 1);
            dynam    = nan(num_pts, num_mols);

            % Initialize counters
            r_steps   = 0;
            rej_steps = 0;

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
                
                % Indexer of reactions in decreading order of propensity
                if mod(r_steps,sort_fr) == 0
                    [~, ind] = sort(prop, "descend");
                end

                % Compute number of allowed reaction firings in time span
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
                            if coef_evol ~= 0
                                r_prop = prop(r);
                                mi_vals(m) = mi_vals(m)+coef_evol*r_prop;
                                si_vals(m) = si_vals(m)+coef_evol^2*r_prop;
                            end
                        end
                    end
                end
                
                % Compute total number of firings in iteration 
                tot_fir = inf; 
                for m=1:num_mols
                    if is_react(m)
                        % Base for both numerators (x^1, x^2 respectively)
                        numtor = max(neg_err*curr_state(m)/(gi_vals(m)),1);
                        % Bases for denominators
                        mi = mi_vals(m);
                        sterm_denom = abs(si_vals(m))-abs(mi^2/tot_prop);
                        % Optimization terms
                        fterm = numtor   / abs(mi);
                        sterm = numtor^2 / sterm_denom;
                        % Minimum for molecule m
                        m_fir = min(fterm,sterm); 
                        % Update if it is the new minimum
                        if m_fir < tot_fir
                            tot_fir = m_fir;
                        end
                    end
                end
                tot_fir = floor(tot_fir * tot_prop); % Num firings in Z+
                
                % Compute alternative num of firings for fast evolving sys
                % Vector with one value per reaction 
                Lj = Inf(1,num_reac); 
                for r = 1:num_reac
                    evol_coefs = mat_evol(r,:);
                    for m = 1:num_mols
                        mol_coef = evol_coefs(m);
                        if mol_coef < 0
                            m_Lj = floor(curr_state(m)/abs(mol_coef));
                            if m_Lj < Lj(r)
                                Lj(r) = m_Lj;
                            end
                        end
                    end
                end
                
                tot_fir_alt = inf;
                for r = 1:num_reac
                    L_alt = (1-theta*(1-(tot_prop/prop(r))))*Lj(r);
                    if L_alt < tot_fir_alt
                        tot_fir_alt = L_alt;
                    end
                end
                tot_fir_alt = floor(tot_fir_alt); % Num firings in Z+

                % Decide num firings among conventional and neg control one
                tot_fir = min(tot_fir_alt, tot_fir);
                tot_fir = max(1,tot_fir);
                
                % Iterate until accepted step (no negative population)
                while true
                    
                    % Condition to break out of loop 
                    accept_step = true;

                    % Sample num of firings per reac from corr binom distr
                    k_fir    = zeros(1,num_reac); 
                    r        = 0;
                    cum_fir  = 0;
                    cum_prop = 0;
    
                    while (cum_fir < tot_fir) && (r < num_reac)
                        r     = r + 1;        % i-th reac in decr prop ord
                        r_ind = ind(r);       % ind of i-th reac in stoich
                        r_prop = prop(r_ind); % propensity of i-th reac
    
                        % Sample num of firings for i-th reac from binom
                        res_fir  = tot_fir  - cum_fir;
                        res_prop = tot_prop - cum_prop;
                        k = binornd(res_fir,r_prop/res_prop);
    
                        cum_fir = cum_fir + k; % Cumulative firings
                        k_fir(r_ind) = k;      % Num firings (stoich ord)
                    end
    
                    % Try evolving the system
                    new_state = curr_state + k_fir * mat_evol;

                    % Check for negative populations 
                    for m=1:num_mols
                        if new_state(m) < 0
                            accept_step = false;
                            rej_steps = rej_steps + 1;
                            tot_fir = max(1,tot_fir/2);
                            break
                        end
                    end
    
                    % Update system
                    if accept_step
                        
                        % Increments R steps counter 
                        r_steps = r_steps + 1;

                        % Compute firing time
                        tau = gamrnd(tot_fir, 1/tot_prop);

                        % Update state and time
                        prev_time  = curr_time;
                        prev_state = curr_state;
                        curr_time  = curr_time + tau;
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

                    end % end sys update

                end % end tot_firings acceptance loop

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
             
            obj.res_s.sim_stat.num_tot_steps = r_steps + rej_steps;
            obj.res_s.sim_stat.num_r_steps = r_steps; 
            obj.res_s.sim_stat.num_rej_steps = rej_steps;
      
        end
        %//////////////////////////////////////////////////////////////////
    
    end

end

classdef SimulationManager < dynamicprops
    % SimulationManager: object for automated execution of simulations
    %     Automatically execute simulations using a grid of parameters. 
    %
    %     Required constructor inputs
    %     grid_file: path to file containing grid of models and algorithms 
    %     verbose: print status update in the console

    properties

        sim_fun_folder     string  = "SimulationFunctions"
        mod_mat_folder     string  = "ModelMatrices"
        par_mat_folder     string  = "ParameterMatrices"
        path_to_grid       string  = "simulation_grid.csv"
        output_folder      string  = "SimulationResults"

    end

    methods

        % RUN SIMULATION GRID /////////////////////////////////////////////
        function run_simulation_grid(obj)
            
            % Attach simulation functions path
            addpath(obj.sim_fun_folder);
            
            % Load simulation grid
            param_grid = readtable(obj.path_to_grid);
            
            % Create general output folder if it does not exist
            if ~exist(obj.output_folder, 'dir')
                mkdir(obj.output_folder)
            end
            
            % Iterate through each line of the input grid
            for line = 1:height(param_grid)
                
                % Retrieve set of objects for simulation
                grid_line = param_grid(line,:);
                
                % Create output folder
                f_path = obj.output_folder+"/"+grid_line.Id(1);
                if ~exist(f_path, 'dir')
                    mkdir(f_path)
                end
                
                % Load model object
                m_path = obj.mod_mat_folder+"/"+grid_line.Model(1)+".mat";
                mod_obj = importdata(m_path);
                
                % Load parameters object
                p_path = obj.par_mat_folder+"/"+grid_line.Params(1)+".mat";
                par_obj = importdata(p_path);
                
                % Create simulation object
                sim_obj = feval(grid_line.Algorithm{1}, mod_obj, par_obj);
                
                % Define result containers sizes
                m_width  = sim_obj.par_s.n_discr_pt;
                m_height = sim_obj.mod_s.num_mols;
                m_depth  = grid_line.Times(1);
                
                % Initialize result containers
                exec_times = nan(m_depth, 1);
                n_tot_step = nan(m_depth, 1);
                discr_time = nan(m_width, m_depth);
                discr_mat  = nan(m_width, m_height, m_depth);

                % Loop for the required amount of simulations left
                for i = 1:m_depth
                    
                    % Run simulation
                    sim_obj.run_simulation()
                    
                    % Store desired results
                    exec_times(i) = sim_obj.res_s.sim_time;
                    n_tot_step(i) = sim_obj.res_s.sim_stat.num_tot_steps;
                    discr_time(:, i) = sim_obj.res_s.full_tim;
                    discr_mat(:,:,i) = sim_obj.res_s.full_mat;

                end

                % Save results
                save(f_path+"/"+"exec_times.mat", "exec_times")
                save(f_path+"/"+"n_tot_step.mat", "n_tot_step")
                save(f_path+"/"+"discr_time.mat", "discr_time")
                save(f_path+"/"+"discr_mat.mat" , "discr_mat" )
                 
            end
            
        end
        % /////////////////////////////////////////////////////////////////
    end
    
end

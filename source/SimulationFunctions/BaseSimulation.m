classdef BaseSimulation < dynamicprops
    % BaseSimulation: parent class for algorithm-specific classes

    properties
        
        mod_s ModStruct
        par_s ParStruct
        res_s ResStruct
        
    end

    methods
        
        % CLASS CONSTRUCTOR ///////////////////////////////////////////////
        function obj = BaseSimulation(mod_s, par_s)

            obj.mod_s = mod_s;
            obj.par_s = par_s;
 
        end
        %//////////////////////////////////////////////////////////////////

        % RUN SIMULATION //////////////////////////////////////////////////
        % Just a placeholder function to overload in the inheriting classes
        function run_algorithm(obj)
            disp("Running " + class(obj) + "(this is just a placeholder)")
        end

        function run_simulation(obj)
            
            % Reset results container
            obj.res_s = ResStruct();
            
            % Start simulation timer
            start_cpu = cputime;

            % Run main algorithm 
            obj.run_algorithm()
            
            % Stop simulation timer and store it
            obj.res_s.sim_time = cputime - start_cpu;
            
        end
        %//////////////////////////////////////////////////////////////////

    end

end

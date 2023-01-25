classdef ModStruct
    
    properties
        
        % Required input fields
        mod_name (1,1) string 
        mol_name (1,:) string
        mat_reag (:,:) double {mustBeInteger, mustBeNonnegative}
        mat_prod (:,:) double {mustBeInteger, mustBeNonnegative}
        st_rates (1,:) double {mustBeNonnegative}
        init_pop (1,:) double {mustBeNonnegative}

        
        % Computed on struct creation
        mat_evol (:,:) double {mustBeInteger}
        num_reac (1,1) double {mustBeInteger}
        num_mols (1,1) double {mustBeInteger}

        % Optionale
        mod_info       string = ''

    end

    methods

        function obj = ModStruct( ...
                mod_name, ...
                mol_name, ...
                mat_reag, ...
                mat_prod, ...
                st_rates, ...
                init_pop)
            
            % Perform size consistency checks
            [num_reactions, num_molecules] = size(mat_reag);

            if length(mol_name) ~= num_molecules
                throw(MException("ModelStructConstructor:Mismatch", ...
                    "Mismatch in number of molecules and names provided"))
            end

            if size(mat_reag) ~= size(mat_prod)
                throw(MException("ModelStructConstructor:Mismatch", ...
                    "Mismatch in sizes of stoichiometric matrices"))
            end

            if length(st_rates) ~= num_reactions
                throw (MException("ModelStructConstructor:Mismatch", ...
                    "Incorrect number of stochastic rates was provided"))
            end 

            if length(init_pop) ~= num_molecules
                throw (MException("BaseSimulatorConstructor:Mismatch", ...
                    "Incorrect num of molecule in initial state vector"))
            end

            % Given from input
            obj.mod_name = mod_name;
            obj.mat_reag = mat_reag;
            obj.mat_prod = mat_prod;
            obj.st_rates = st_rates;
            obj.init_pop = init_pop;
            
            % Computed from input
            obj.mat_evol = mat_prod - mat_reag;
            obj.num_reac = num_reactions;
            obj.num_mols = num_molecules;

        end

    end

end

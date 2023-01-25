classdef ParStruct
    
    properties
        
        % Required input fields
        t_max        (1,1) double {mustBeNonnegative}
        n_discr_pt   (1,1) double {mustBeNonnegative}

        % Algorithm specific fields
        n_crit_thr   (1,1) double {mustBeNonnegative}
        epsilon      (1,1) double {mustBeNonnegative}
        n_ssa_steps  (1,1) double {mustBeNonnegative}
        reord_freq   (1,1) double {mustBeNonnegative}
        theta        (1,1) double {mustBeNonnegative}

    end

    methods

        function obj = ParStruct(t_max)

            obj.t_max = t_max;

        end

    end

end

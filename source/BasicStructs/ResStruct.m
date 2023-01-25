classdef ResStruct

    properties
       
        sim_time (1,1) double {mustBeNonnegative}
        full_mat (:,:) double {mustBeNonnegative}
        full_tim (:,:) double {mustBeNonnegative}
        disc_mat (:,:) double {mustBeNonnegative}

        sim_stat = struct( ...
            "num_tot_steps", 0, ...
            "num_ssa_steps", 0, ...
            "num_tau_steps", 0, ...
            "num_r_steps"  , 0, ...
            "num_s_steps"  , 0, ...
            "num_rej_steps", 0  )

    end

end
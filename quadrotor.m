%% QUADROTOR quadrotor simulator class %
%
% 
%
classdef quadrotor
    
    %% quadrotor public variables
    properties (Constant, GetAccess='public')

    end
    
    %% quadrotor private variables
    properties (GetAccess='private', SetAccess='private')
        %% propeller properties
        % propeller moment constant vector (this should be 4x1)
        km = 2.400e-7*[1;-1;1;-1];
        % propeller force constant (usually negative 4x1)
        kf = 9.169e-6*[-1; -1; -1; -1];
        % propeller inertial matrix
        Jpz = 3.46e-3*[1; 1; 1; 1];
        Jpx = 0*[1; 1; 1; 1];
        %% body inertial/geometry properties
        % total mass
        m = 0.538;
        % fuselage inertia matrix
        J = [0.0084728 0 0; 0 0.00845974 0; 0 0 0.00811587];
        % semi-side of quadrotor in meters (quad has area 4L^2)
        L = 0.0825;
        %% state variables
        % position in NED coordinates (m)
        p = zeros(3,1);
        % velocity wrt to local in NED coordinates (m/s)
        v = zeros(3,1);
        % angular velocity wrt to local in body coordinates (rad/s)
        w = zeros(3,1);
        % quaternion of fuselage wrt local
        q = [1; 0; 0; 0];
        %% physical constants
        g = 9.81;
        
    end
    
    %% Navigator public methods
    methods (Access=public)
        
        function obj = update_state(obj,u,dt)
            %% UPDATE_STATE: updates system state based on new inputs
            %
            % inputs:  u  - 4x1 Real - propeller speeds with respect to fuselage (rad/s)
            %          dt - 1x1 Real - sample period (secs)
            % outputs: none
            %
            
            x = [obj.p; obj.v; obj.w; obj.q];
            dx_dt = obj.state_derivative(u);
            x = x + dx_dt * dt;
            obj.p = x(1:3);
            obj.v = x(4:6);
            obj.w = x(7:9);
            obj.q = x(10:13);
            obj.q = obj.q/norm(obj.q,2); % quaternion normalization here
            
        end
        
        function p = get_position_ned(obj)
            p = obj.p;
        end
        
        function v_ned = get_velocity_ned(obj)
            v_ned = obj.v;
        end
        
        function w = get_ang_vel_body(obj)
            w = obj.w;
        end
        
        function q = get_quaternion(obj)
            q = obj.q;
        end
        
    end
    
    %% quadrotor private methods
    methods (Access=private)
        
        function x_dot = state_derivative(obj,u)
            %% STATE_DERIVATIVE: computes state derivative x_dot = f(x,u)
            %
            % inputs:  u  - 4x1 Real - propeller speeds with respect to fuselage (rad/s)
            % outputs: x_dot - 13x1 Real - derivative of system state
            %
            % notice: x = (p v w q)
            
            p_dot = obj.v;
            
            D = quat2dcm(obj.q');
            v_dot = 1/obj.m*D'*obj.fuselage_force(u);
            
            M = obj.fuselage_moment(u);
            w_dot = (obj.J)\(M - obj.cross_operator(obj.w)*obj.J*obj.w);
            
            A = [0 -obj.w'; obj.w -obj.cross_operator(obj.w) ];
            q_dot = 1/2*A*obj.q;
            
            x_dot = [ p_dot; v_dot; w_dot; q_dot ];
            
        end
        
        function Mi = prop_internal_moment(obj, u, id)
            %% PROP_INTERNAL_MOMENT: computes propeller internal (prop-fuselage) moment
            %
            % inputs:  u  - 4x1 Real - propeller speeds with respect to fuselage (rad/s)
            %          id - 1x1 Real - propeller id number in {1, 2, 3, 4}
            % outputs: Mi - 3x1 Real - interaction moment in fuselage coordinates (Nm)
            %
            % notice: wi is positive when pointing in the direction of Z in
            % fuselage reference frame!!!
            P = obj.w(1); Q = obj.w(2); R = obj.w(3);
            wi = u(id);
            Mi = [0; 0; obj.km(id)*wi^2] + (obj.Jpz(id)-obj.Jpx(id))*(R+wi)*[-Q; P; 0];
        end
        
        function Fi = prop_internal_force(obj, u, id)
            %% PROP_INTERNAL_FORCE: computes propeller internal (prop-fuselage) force
            %
            % inputs:  u  - 4x1 Real - propeller speeds with respect to fuselage (rad/s)
            %          id - 1x1 Real - propeller id number in {1, 2, 3, 4}
            % outputs: Fi - 3x1 Real - interaction force in fuselage coordinates (N)
            %
            % notice: kf is such that it is normally negative!!
            wi = u(id);
            Fi = [0; 0; obj.kf(id)*wi^2];
        end
        
        function Fb = fuselage_force(obj, u)
            %% FUSELAGE_FORCE: computes fuselage resultant force
            %
            % inputs:  u  - 4x1 Real - propeller speeds with respect to fuselage (rad/s)
            % outputs: Fb - 3x1 Real - fuselage resultant force in body coordinates (N)
            %
            
            % computes dcm from local to body
            D = quat2dcm(obj.q');
            % weight force in NED coordinates
            Wl = obj.m*obj.g*[0;0;1];
            % resultant force computation
            Fb = D*Wl;
            for id = 1:4
                Fb = Fb + obj.prop_internal_force(u, id);
            end
            
        end
        
        function Mb = fuselage_moment(obj, u)
            %% FUSELAGE_MOMENT: computes fuselage resultant moment
            %
            % inputs:  u  - 4x1 Real - propeller speeds with respect to fuselage (rad/s)
            % outputs: Mb - 3x1 Real - fuselage resultant moment in body coordinates (Nm)
            %
            
            Mb = zeros(3,1);
            for id=1:4
                Mb = Mb + obj.prop_internal_moment(u, id);
                Fi = obj.prop_internal_force(u, id);
                Pi = obj.pos_prop_wrt_cg(id);
                Mb = Mb + cross(Pi,Fi);
            end
            
            
            
        end
        
        function p = pos_prop_wrt_cg(obj, id)
            %% POS_PROP_WRT_CG: computes position of prop_i wrt to center of mass
            %
            % outputs: p  - 3x1 Real - position of prop_i wrt to center of mass (m)
            %
            
            % semi-side of quadrotor
            l = obj.L;
            
            switch id
                case 1
                    p = [l;l;0];
                case 2
                    p = [-l;l;0];
                case 3
                    p = [-l;-l;0];
                case 4
                    p = [l;-l;0];
                otherwise
                    error('pos_prop_wrt_cg(): please choose a valid propeller id number!!');
            end
            
        end
        
        function Vx = cross_operator(~, v)
            P = v(1); Q = v(2); R = v(3);
            Vx = [0 -R Q; R 0 -P; -Q P 0];
        end
        
    end
    
    methods (Static)
        
        function [thrust_sp, att_control, rates_sp_prev_out, rates_prev_out, rates_int_out] = pixhawk_mc_att_control(w_m, q_r, q_m, rates_sp_prev, rates_prev, rates_int, dt)
            
            rates_int_out = zeros(3,1);
            
            %% pixhawk control configuration parameters
            
            % common gain
            Kc = 3;
            % attitude proportional gains (inner loop)
            ang_P_gains = Kc*[5; 4; 3];
            % attitude proportional gains (outter loop)
            rate_P_gains = Kc*[0.15; 0.18; 0.40];
            rate_I_gains = Kc*[0; 0; 0];
            rate_D_gains = Kc*[0; 0; 0];
            % integral control limit
            RATES_I_LIMIT = 100000000;
            % yaw feedforward weight
            yaw_ff = 0;
            % rate feedforward
            rate_ff = [0;0;0];
            % yaw rate setpoint from pilot
            yaw_sp_move_rate = 0;
            % minimum takeoff thrust
            MIN_TAKEOFF_THRUST = sqrt(0.538*9.81/4/9.169e-6);
            
            
            %% pixhawk inner attitude control law (mc_att_control_main.cpp)
            
            thrust_sp = sqrt(0.538*9.81/4/9.169e-6); % pilot value here!!! now we are only interested in attitude and will disregard this
            
            % construct attitude setpoint rotation matrix */
            R_sp = quat2dcm(q_r')';
            % get current rotation matrix from control state quaternions */
            R = quat2dcm(q_m')';
            
            % all input data is ready, run controller itself */
            
            % try to move thrust vector shortest way, because yaw response is slower than roll/pitch */
            R_z = R * [0; 0; 1];
            R_sp_z = R_sp * [0; 0; 1];
            
            % axis and sin(angle) of desired rotation */
            e_R = R' * cross(R_z, R_sp_z);
            
            % calculate angle error */
            e_R_z_sin = norm(e_R,2);
            e_R_z_cos = dot(R_z, R_sp_z);
            
            % calculate weight for yaw control */
            yaw_w = R_sp(3,3) * R_sp(3,3);
            
            % calculate rotation matrix after roll/pitch only rotation */
            
            if e_R_z_sin > 0
                % get axis-angle representation */
                e_R_z_angle = atan2(e_R_z_sin, e_R_z_cos);
                e_R_z_axis = e_R / e_R_z_sin;
                
                e_R = e_R_z_axis * e_R_z_angle;
                
                % cross product matrix for e_R_axis */
                e_R_cp = zeros(3,3);
                e_R_cp(1,2) = -e_R_z_axis(3);
                e_R_cp(1,3) = e_R_z_axis(2);
                e_R_cp(2,1) = e_R_z_axis(3);
                e_R_cp(2,3) = -e_R_z_axis(1);
                e_R_cp(3,1) = -e_R_z_axis(2);
                e_R_cp(3,2) = e_R_z_axis(1);
                
                % rotation matrix for roll/pitch only rotation */
                % this was weird pixhawk:
                %R_rp = R * (eye(3) + e_R_cp * e_R_z_sin + e_R_cp * e_R_cp * (1 - e_R_z_cos));
                R_rp = R * (e_R_z_cos * eye(3) + e_R_cp * e_R_z_sin + e_R_z_axis * e_R_z_axis' * (1 - e_R_z_cos));
            else
                % zero roll/pitch rotation */
                R_rp = R;
            end
            
            % R_rp and R_sp has the same Z axis, calculate yaw error */
            R_sp_x = R_sp*[1;0;0];
            R_rp_x = R_rp*[1;0;0];
            % this was weird pixhawk:
            %e_R(3) = atan2(dot(cross(R_rp_x,R_sp_x),R_sp_z), dot(R_rp_x,R_sp_x)) * yaw_w;
            e_R(3) = atan2(dot(cross(R_rp_x,R_sp_x),R_sp_z), dot(R_rp_x,R_sp_x));
            
            if e_R_z_cos < 0
                % for large thrust vector rotations use another rotation method:
                %  calculate angle and axis for R -> R_sp rotation directly */
                q = dcm2quat(R_sp'*R)';
                e_R_d = q(2:4);
                e_R_d = e_R_d/norm(e_R_d,2);
                e_R_d = e_R_d * 2 * atan2(norm(e_R_d,2), q(1));
                
                % use fusion of Z axis based rotation and direct rotation */
                direct_w = e_R_z_cos * e_R_z_cos * yaw_w;
                e_R = e_R * (1 - direct_w) + e_R_d * direct_w;
            end
            
            % calculate angular rates setpoint */
            rates_sp = ang_P_gains .* e_R;
            
            %	/* limit rates */ THIS PART HAS NOT BEEN TRANSLATED FROM C CODE, MAYBE LATER?
            %	for (int i = 0; i < 3; i++) {
            %		if (_v_control_mode.flag_control_velocity_enabled && !_v_control_mode.flag_control_manual_enabled) {
            %			_rates_sp(i) = math::constrain(_rates_sp(i), -_params.auto_rate_max(i), _params.auto_rate_max(i));
            %		} else {
            %			_rates_sp(i) = math::constrain(_rates_sp(i), -_params.mc_rate_max(i), _params.mc_rate_max(i));
            %		}
            %	}
            
            % feed forward yaw setpoint rate */
            rates_sp(3) = rates_sp(3) + yaw_sp_move_rate * yaw_w * yaw_ff;
            
            %% pixhawk outer attitude control law (mc_att_control_main.cpp)
            
            % current body angular rates */
            rates = w_m;
            
            % angular rates error */
            rates_err = rates_sp - rates;
            att_control = rate_P_gains.*rates_err + rate_D_gains.*(rates_prev-rates)/dt + rates_int + rate_ff.*(rates_sp-rates_sp_prev)/dt;
            rates_sp_prev_out = rates_sp;
            rates_prev_out = rates;
            
            % update integral only if not saturated on low limit and if motor commands are not saturated */
            if (thrust_sp > MIN_TAKEOFF_THRUST) % i suppressed the saturation here!! we can put it back later
                for i=1:3
                    if (abs(att_control(i)) < thrust_sp)
                        rate_i = rates_int(i) + rate_I_gains(i) * rates_err(i) * dt;
                        
                        if (rate_i > -RATES_I_LIMIT && rate_i < RATES_I_LIMIT && att_control(i) > -RATES_I_LIMIT && att_control(i) < RATES_I_LIMIT)
                            rates_int_out(i) = rate_i;
                        end
                    end
                end
            end
            
        end
        
    end
    
end




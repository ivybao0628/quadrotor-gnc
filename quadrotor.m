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
            obj.p = 0*x(1:3);
            obj.v = 0*x(4:6);
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
    
    %% Navigator private methods
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
    
end




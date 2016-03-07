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
        R_rp = R * (eye(3) + e_R_cp * e_R_z_sin + e_R_cp * e_R_cp * (1 - e_R_z_cos));
	else
		% zero roll/pitch rotation */
		R_rp = R;
	end

	% R_rp and R_sp has the same Z axis, calculate yaw error */
    R_sp_x = R_sp*[1;0;0];
    R_rp_x = R_rp*[1;0;0];
	e_R(3) = atan2(dot(cross(R_rp_x,R_sp_x),R_sp_z), dot(R_rp_x,R_sp_x)) * yaw_w;

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
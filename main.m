%% main simulation file
%
%
%

clear all;
close all;

%% drone init
drone = quadrotor();

%% time parameters
dt = 0.01;
t0 = 0;
tF = 60;
t = t0:dt:tF;

%% animation configuration
anim_env = Aero.Animation;
anim_env.FramesPerSecond = 12;
anim_env.TimeScaling = 1;
dt_frame = 1/anim_env.FramesPerSecond;
dt_elapsed = 0;
n_frames = 0;
anim_body = anim_env.createBody('pa24-250_orange.ac','Ac3d');
anim_env.Bodies{1}.TimeSeriesSourceType = 'Array6DoF';
anim_env.Bodies{1}.TimeSeriesSource = zeros(length(t),7);


%% pixhawk init
rates_int = zeros(3,1);
rates_sp_prev = zeros(3,1);
rates_prev = zeros(3,1);

%% logger init 
p_log = zeros(3,length(t));
v_log = zeros(3,length(t));
w_log = zeros(3,length(t));
q_log = zeros(4,length(t));
psi_r_log = zeros(1,length(t));
the_r_log = zeros(1,length(t));
phi_r_log = zeros(1,length(t));
psi_m_log = zeros(1,length(t));
the_m_log = zeros(1,length(t));
phi_m_log = zeros(1,length(t));
att_control_log = zeros(3,length(t));

%% main loop
for ii=1:length(t)
    
    %% measured values
    p_m = drone.get_position_ned();
    v_m = drone.get_velocity_ned();
    w_m = drone.get_ang_vel_body();
    q_m = drone.get_quaternion();
    [psi_m, the_m, phi_m] = quat2angle(q_m');
    
    %% desired values
    p_r = zeros(3,1);
    v_r = zeros(3,1);
    w_r = zeros(3,1);
    psi_r = 0*pi/180+0*2*pi/tF*t(ii);
    the_r = 0*pi/180+0*2*pi/tF*t(ii);
    phi_r = 15*pi/180+0*2*pi/tF*t(ii);
    q_r = angle2quat( psi_r, the_r, phi_r )';
    
    %% pixhawk control law
    [thrust_sp, att_control, rates_sp_prev, rates_prev, rates_int_out] = mod_pixhawk_mc_att_control(w_m, q_r, q_m, rates_sp_prev, rates_prev, rates_int, dt);
    
    %% update drone state according to control law and dyamics model
    u0 = [-1;+1;-1;+1] * thrust_sp;
    uP = [+1;-1;-1;+1] * att_control(1);
    uQ = [-1;-1;+1;+1] * att_control(2);
    uR = [-1;-1;-1;-1] * att_control(3);
    
    drone = drone.update_state(u0+uP+uQ+uR,dt);
    
    %% save data for plotting later
    p_log(:,ii) = drone.get_position_ned();
    v_log(:,ii) = drone.get_velocity_ned();
    w_log(:,ii) = drone.get_ang_vel_body();
    q_log(:,ii) = drone.get_quaternion();
    psi_r_log(1,ii) = psi_r;
    the_r_log(1,ii) = the_r;
    phi_r_log(1,ii) = phi_r;
    psi_m_log(1,ii) = psi_m;
    the_m_log(1,ii) = the_m;
    phi_m_log(1,ii) = phi_m;
    att_control_log(:,ii) = att_control;
    
    %% animation constructions
    dt_elapsed = dt_elapsed + dt;
    if ( dt_elapsed > dt_frame ) 
        dt_elapsed = 0;
        n_frames = n_frames + 1;
        anim_env.Bodies{1}.TimeSeriesSource(n_frames,1) = t(ii);
        [psi, the, phi] = quat2angle(drone.get_quaternion()');
        anim_env.Bodies{1}.TimeSeriesSource(n_frames,5:7) = [phi the psi]; 
    end
    
    
end

%% fix and plot animation
anim_env.Bodies{1}.TimeSeriesSource = anim_env.Bodies{1}.TimeSeriesSource(1:n_frames,:);
anim_env.show();
anim_env.play()

%% plot results
figure;
subplot(3,1,1);
plot(t, 180/pi*psi_r_log, t, 180/pi*psi_m_log);
ylabel('Yaw (deg)');
subplot(3,1,2);
plot(t, 180/pi*the_r_log, t, 180/pi*the_m_log);
ylabel('Pitch (deg)');
subplot(3,1,3);
plot(t, 180/pi*phi_r_log, t, 180/pi*phi_m_log);
ylabel('Roll (deg)');
xlabel('time (sec)');








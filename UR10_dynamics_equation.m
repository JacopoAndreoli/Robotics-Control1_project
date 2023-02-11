
% Jacopo Andreoli n° matrc: 2011655

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UR10 Derivation of the kinematics and dynamics equations  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
%UR10 Kinematics %
%%%%%%%%%%%%%%%%%%



%% import casadi and define robot parameters
import casadi.*

theta_list = [0, 0, 0, 0, 0, 0];
d_list = [0.1273, 0, 0, 0.163941, 0.1157, 0.0922];
alpha_list = [pi/2, 0, 0, pi/2, -pi/2, 0];
a_list = [0, -0.612, -0.5723, 0, 0, 0];

% the dynamic parameters associated to the first link are set to zero
% since the link does not influence the dynamics of the manipulator

% define the links center of mass coordinates
p_matrix = [ 0, 0, 0;
            0.021,  0.0, 0.027;
            0.38, 0.0, 0.158;
            0.24, 0.0, 0.068;
            0.0,  0.007, 0.018;
            -0.0, 0.007, 0.018;
            0.0, 0.0, -0.026]; 

% define a vector with the joint type (if true the joint is revolute)
num_dof = 6;
joint_type_list = [true, true, true, true, true, true];



%% symbolic variables

% joint positions
q = SX.sym('q', [num_dof,1]);
% joint coordinates
q_dot = SX.sym('q_dot', [num_dof,1]);



%% compute the expressions of the DH reference frames

% allocate the space
R_DH = cell(num_dof,1);
l_DH = cell(num_dof,1);

%iterate all the joints
for j = 1 : num_dof
    
    % check the joint type and determine d_j and theta_j
    if joint_type_list(j)
        % if the joint is rev theta_j = theta_j^0 + q_j
        theta_j = theta_list(j) + q(j);
        d_j = d_list(j);
    else
        % if the joint is prismatic d_j = d_j^0 + q_j
        theta_j = theta_list(j);
        d_j = d_list(j) + q(j);
    end
    
    % get DH transformation
    [R_DH{j}, l_DH{j}] = DH_transform(alpha_list(j), a_list(j), d_j, theta_j);
    
end

% allocate variables
R_j_0_list = cell(num_dof+1,1);
l_j_0_list = cell(num_dof+1,1);

% initialize the variables
l_j_0_list{1} = SX.zeros(3,1);
R_j_0_list{1} = SX.eye(3);

%compute the relative transformations
for j = 2 : num_dof+1
   % l_j_0 = l_j-1_0 + R_j-1^0*l_j^j-1
   l_j_0_list{j} = l_j_0_list{j-1} + R_j_0_list{j-1}*l_DH{j-1};
   % R_j_0 = R_j-1^0*R_j^j-1
   R_j_0_list{j} = R_j_0_list{j-1}*R_DH{j-1};
end

%% convert R_e_0 in yaw pitch and roll angles

[yaw, pitch, roll] = rot_2_YPR_angles(R_j_0_list{num_dof+1});

%% get the numeric function and the mex function

% get functions evaluable numerically
l_matrix = [l_j_0_list{2},...
            l_j_0_list{3},...
            l_j_0_list{4},...
            l_j_0_list{5},...
            l_j_0_list{6},...
            l_j_0_list{7}];
R_matrix = [R_j_0_list{2},...
            R_j_0_list{3},...
            R_j_0_list{4},...
            R_j_0_list{5},...
            R_j_0_list{6},...
            R_j_0_list{7}];
        
f_kin_l = Function('f_kin_l',{q}, {l_matrix});
f_kin_R = Function('f_kin_R',{q}, {R_matrix});
f_x = Function('f_x',{q},{[l_j_0_list{num_dof+1}; yaw; pitch; roll]});

% save the casadi function
save('f_kin_l','f_kin_l')
save('f_kin_R','f_kin_R')
save('f_x','f_x')

% generate the mex function
opts = struct('main', true,...
              'mex', true);
f_kin_l.generate('f_kin_l_mex.c',opts);
f_kin_R.generate('f_kin_R_mex.c',opts);
f_x.generate('f_x_mex.c',opts);
mex f_kin_l_mex.c
mex f_kin_R_mex.c
mex f_x_mex.c


%% compute the positions of the centers of mass

% allocate variables
P_p = cell(num_dof+1,1);

% compute positions
for j = 1 : num_dof+1
    P_p{j} = l_j_0_list{j} + R_j_0_list{j}*p_matrix(j,:)';
end

%% express axis orientation through yaw,pitch,roll angles

% allocate variables
yaw_angles = cell(num_dof+1,1);
pitch_angles = cell(num_dof+1,1);
roll_angles = cell(num_dof+1,1);
alpha_vect{1} = 0;

% compute Z angle
for j = 2 : num_dof+1
    [yaw_angles{j}, pitch_angles{j}, roll_angles{j}] = rot_2_YPR_angles(R_j_0_list{j});
end

%here we define the matrix that express relation between analytical and geometric
%jacobian; the metrix is equivalent to the rotation matrix needed to
%express the angular velocity given the Euler angles velocity
%rappresentation

T_phi = cell(3,3);

for j = 2 : num_dof+1
   T_phi{j} = [[1, 0, sin(pitch_angles{j})];...
               [0, cos(roll_angles{j}), -sin(roll_angles{j})*cos(pitch_angles{j})];...
               [0, sin(roll_angles{j}), cos(roll_angles{j})*cos(pitch_angles{j})]];
end

%% get jacobians


%here we exploit the computation of the analitical Jacobian defined through
%the Casadi function and the we translate it into its associated geometric
%Jacobian through the Matrix T_phi evaluated above

% allocate the space
J_p = cell(num_dof+1,1);
J_o = cell(num_dof+1,1);

for j = 2 : num_dof+1
    % linear velocities
    J_p{j} = jacobian(P_p{j}, q); % for linear velocity, geometric and analitycal Jacobian coincide
    
    % angular velocties
    J_o{j} = T_phi{j}*jacobian([yaw_angles{j}, pitch_angles{j}, roll_angles{j}], q);
end



%% compute the inertia matrix and the potential energy

% init the matrix with the first link - the tensor is always referred to
% the link reference frame; this allow us to consider it a constant
% quantity
B = m_vect(2)*J_p{2}'*J_p{2} + J_o{2}'*R_j_0_list{2}'*I_tensor{2}*R_j_0_list{2}*J_o{2};
g = [0;0;-9.80665];
U = -m_vect(2)*g'*P_p{2};

for j = 3 : num_dof+1
    B = B + m_vect(j)*J_p{j}'*J_p{j} + J_o{j}'*R_j_0_list{j}'*I_tensor{j}*R_j_0_list{j}*J_o{j};
    U = U -m_vect(j)*g'*P_p{j};
end 

%% compute dynamics equations

% terms dependent on velocities
H = cell(num_dof,1);
for i=1:num_dof
    % get H1 and H2
    H{i} = jacobian(B(i,:),q) - 0.5*reshape(jacobian(B,q(i)), num_dof, num_dof);
end
C_vect = [q_dot'*H{1}*q_dot;
          q_dot'*H{2}*q_dot;
          q_dot'*H{3}*q_dot;
          q_dot'*H{4}*q_dot;
          q_dot'*H{5}*q_dot;
          q_dot'*H{6}*q_dot;];

% gravitational contribution
g_vect = jacobian(U, q);

%% get functions evaluable numerically

% get casadi functions
f_B = Function('f_B',{q}, {B});
f_U = Function('f_U',{q}, {U});
f_C_vect = Function('f_C_vect',{q, q_dot}, {C_vect});
f_g_vect = Function('f_g_vect',{q}, {g_vect});

% get mex functions
opts = struct('main', true,...
              'mex', true);
f_B.generate('f_B_mex.c',opts);
mex f_B_mex.c
f_U.generate('f_U_mex.c',opts);
mex f_U_mex.c
f_C_vect.generate('f_C_vect_mex.c',opts);
mex f_C_vect_mex.c
f_g_vect.generate('f_g_vect_mex.c',opts);
mex f_g_vect_mex.c


%% evaluation of POS

POS = cell(100,3);
ANG = cell(100,3);
for i=1:100
    x = f_x_mex('f_x', Q_kin_test(i,:));
    POS{i,1} = x(1);
    POS{i,2} = x(2);
    POS{i,3} = x(3);
    ANG{i,1} = x(4);
    ANG{i,2} = x(5);
    ANG{i,3} = x(6);
end


TAU = cell(100,6);
for i=1:1000
    B_i_q_ddot = full(f_B_mex('f_B',Q_dyn_test(i,:)))*(Q_ddot_dyn_test(i,:)');
    C_dot = full(f_C_vect_mex('f_C_vect', Q_dyn_test(i,:), Q_dot_dyn_test(i,:)));
    G_q = -full(f_g_vect_mex('f_g_vect', Q_dyn_test(i,:)))';
    tau = (B_i_q_ddot + C_dot + G_q)';
    TAU{i,1} = tau(1);
    TAU{i,2} = tau(2);
    TAU{i,3} = tau(3);
    TAU{i,4} = tau(4);
    TAU{i,5} = tau(5);
    TAU{i,6} = tau(6);
    
end


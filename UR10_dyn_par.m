%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UR10 Derivation of the dynamics equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% import casadi and define robot parameters
import casadi.*

% links mass
m_vect = [0, 7.1, 12.7, 4.27, 2, 2, 0.365];

% links length and radius
l_vect = [0, 0.1273, 0.612, 0.5723, 0.163941, 0.1157, 0.0922];
r_vect = [0, 0.2, 0.15, 0.15, 0.1, 0.1, 0.1];
         
% define the links inertia tensor
% (expressed in the correspondent DH reference frame)
I_tensor = cell(7,1);
I_tensor{1} = diag([0,0,0]);
% link 1 (y parallel to the cylender axis)
I_tensor{2} = diag([m_vect(1)/12*(l_vect(1)^2+3*r_vect(1)^2),...
                    0.5*m_vect(1)*r_vect(1)^2,...
                    m_vect(1)/12*(l_vect(1)^2+3*r_vect(1)^2)]);
% link 2 (x parallel to the cylender axis)
I_tensor{3} = diag([0.5*m_vect(2)*r_vect(2)^2,...
                    m_vect(2)/12*(l_vect(2)^2+3*r_vect(2)^2),...
                    m_vect(2)/12*(l_vect(2)^2+3*r_vect(2)^2)]);
% link 3 (x parallel to the cylender axis)
I_tensor{4} = diag([0.5*m_vect(3)*r_vect(3)^2,...
                    m_vect(3)/12*(l_vect(3)^2+3*r_vect(3)^2),...
                    m_vect(3)/12*(l_vect(3)^2+3*r_vect(3)^2)]);
% link 4 (y parallel to the cylender axis)
I_tensor{5} = diag([m_vect(4)/12*(l_vect(4)^2+3*r_vect(4)^2),...
                    0.5*m_vect(4)*r_vect(4)^2,...
                    m_vect(4)/12*(l_vect(4)^2+3*r_vect(4)^2)]);
% link 5 (y parallel to the cylender axis)
I_tensor{6} = diag([m_vect(5)/12*(l_vect(5)^2+3*r_vect(5)^2),...
                    0.5*m_vect(5)*r_vect(5)^2,...
                    m_vect(5)/12*(l_vect(5)^2+3*r_vect(5)^2)]);
% link 6 (z parallel to the cylender axis)
I_tensor{7} = diag([m_vect(6)/12*(l_vect(6)^2+3*r_vect(6)^2),...
                    m_vect(6)/12*(l_vect(6)^2+3*r_vect(6)^2),...
                    0.5*m_vect(6)*r_vect(6)^2]);
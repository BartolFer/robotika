#!/usr/bin/env python

TESTING = False;

import sys
import datetime
from math import pi, floor, tau
from math import sin, cos
from math import sqrt
from math import atan2, asin, acos
import numpy as np
import copy
if not TESTING:
    import csv
    import rospy
    import geometry_msgs.msg
    from std_msgs.msg import String
    from std_msgs.msg import Header
    from geometry_msgs.msg import Pose
    from geometry_msgs.msg import PoseStamped
    from geometry_msgs.msg import Quaternion
    from geometry_msgs.msg import PoseArray
    from moveit_msgs.msg import RobotState
    from sensor_msgs.msg import JointState
    from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
    import tf
    from open_manipulator_msgs.srv import SetJointPosition, SetJointPositionRequest
pass
from numpy import genfromtxt

#import pandas as pd

class Polynome(list):
    def __call__(self, t):
        x = 1;
        y = 0;
        for a in self:
            y += a * x;
            x *= t;
        pass
        return y;
    pass
    def derivative(self): return Polynome(i * self[i] for i in range(1, len(self)));
pass
class LinkedNode:
    _last_node = None;
    prev = next = None;
    def __init__(self, data):
        self.data = data;
        LinkedNode._last_node = self;
    pass
    def insertBefore(self, data):
        node = LinkedNode(data);
        if self.prev:
            self.prev.next = node;
        pass
        node.prev = self.prev;
        node.next = self;
        self.prev = node;
        return node;
    pass
    def insertAfter(self, data):
        node = LinkedNode(data);
        if self.next:
            self.next.prev = node;
        pass
        node.prev = self;
        node.next = self.next;
        self.next = node;
        return node;
    pass

    def _count(self):
        while self.prev: self = self.prev;
        c = 1;
        while self.next:
            c += 1;
            self = self.next;
        pass
        return c;
    pass
    def _toList(self):
        while self.prev: self = self.prev;
        l = [];
        while self.next:
            l.append(self.data);
            self = self.next;
        pass
        l.append(self.data);
        return l;
    pass
pass

from functools import wraps;
@wraps(print)
def tprint(*args, **kw):
    if TESTING: print(*args, **kw);
pass

class RoboMath:
    def __init__(self):
        # Variables
        self._q = [0, 0, 0, 0]

        # Parameters
        self.l = [0.077, 0.128, 0.024, 0.124, 0.044, 0.105]
        self.gamma = atan2(self.l[1], self.l[2])
        self.h = (self.l[1]**2 + self.l[2]**2)**0.5
    pass

    # LAB 1
    def get_dk(self, q):
        # Implement here direct kinematics
        # INPUT: q as a vector 4x1
        # OUTPUT: w as a vector 6x1

        # print ('Enter DK - q: ', q)

        self.DH_params =  [[q[0],               self.l[0],  0,          -pi/2], 
                           [q[1]-self.gamma,    0,          self.h,     0],
                           [q[2]+self.gamma,    0,          self.l[3],  0],
                           [q[3],               0,          self.l[4],  -pi/2],
                           [0,                  self.l[5],  0,          0]]

        T_t = np.eye(4)

        for i in range(0, 5):
            T_temp = np.eye(4)
            T_temp[0][0] =  cos(self.DH_params[i][0])
            T_temp[0][1] = -cos(self.DH_params[i][3])*sin(self.DH_params[i][0])
            T_temp[0][2] =  sin(self.DH_params[i][3])*sin(self.DH_params[i][0])
            T_temp[0][3] =  self.DH_params[i][2]*cos(self.DH_params[i][0])
            T_temp[1][0] =  sin(self.DH_params[i][0])
            T_temp[1][1] =  cos(self.DH_params[i][3])*cos(self.DH_params[i][0])
            T_temp[1][2] = -sin(self.DH_params[i][3])*cos(self.DH_params[i][0])
            T_temp[1][3] =  self.DH_params[i][2]*sin(self.DH_params[i][0])
            T_temp[2][1] =  sin(self.DH_params[i][3])
            T_temp[2][2] =  cos(self.DH_params[i][3])
            T_temp[2][3] =  self.DH_params[i][1]

            #print (T_temp)

            T_t = copy.deepcopy(np.matmul(T_t, T_temp))

        w = np.zeros(6)
        w[0] = T_t[0][3]
        w[1] = T_t[1][3]
        w[2] = T_t[2][3]
        w[3] = T_t[0][2]
        w[4] = T_t[1][2]
        w[5] = T_t[2][2]

        return w

    def get_ik(self, w):
        # Implement here inverse kinematics
        # INPUT (1): w 6x1 as a tool configuration vector, w = [x, y, z, wx, wy, wz]
        # INPUT (2): q0 4x1 as a joint_state vector of temp robot's position
        # OUTPU: q 4x1 as a closest feasible inverse solution to q0

        # Get q1
        s_q1 = atan2(w[1] , w[0] )

        # Get q3
        star_1 = w[0]*cos(s_q1) + w[1]*sin(s_q1)     - self.l[5]*(w[3]*cos(s_q1) + w[4]*sin(s_q1)) + self.l[4]*w[5]
        star_2 = -(w[2] - self.l[0] - self.l[5]*w[5] - self.l[4]*(w[3]*cos(s_q1) + w[4]*sin(s_q1)) )
        help_q3 = acos( (1/(2*self.l[3]*self.h) )*( ( star_1 )**2 + ( star_2 )**2 - self.l[3]**2 - self.h**2))     

        s_q3_1 = -self.gamma + help_q3
        s_q3_2 = -self.gamma - help_q3
        
        # Get q2
        help_q2_a1 = star_1
        help_q2_b1_1 = self.l[3]*cos(s_q3_1) + self.h*cos(self.gamma)
        help_q2_c1_1 = self.l[3]*sin(s_q3_1) - self.h*sin(self.gamma)
        help_q2_b1_2 = self.l[3]*cos(s_q3_2) + self.h*cos(self.gamma)
        help_q2_c1_2 = self.l[3]*sin(s_q3_2) - self.h*sin(self.gamma)
        help_q2_a2 = star_2
        help_q2_b2_1 = self.l[3]*sin(s_q3_1) - self.h*sin(self.gamma)
        help_q2_c2_1 = self.l[3]*cos(s_q3_1) + self.h*cos(self.gamma)
        help_q2_b2_2 = self.l[3]*sin(s_q3_2) - self.h*sin(self.gamma)
        help_q2_c2_2 = self.l[3]*cos(s_q3_2) + self.h*cos(self.gamma)      

        s_q2_1 = atan2(help_q2_a2*help_q2_b1_1 - help_q2_a1*help_q2_b2_1, help_q2_a2*help_q2_c1_1 + help_q2_a1*help_q2_c2_1)
        s_q2_2 = atan2(help_q2_a2*help_q2_b1_2 - help_q2_a1*help_q2_b2_2, help_q2_a2*help_q2_c1_2 + help_q2_a1*help_q2_c2_2)

        # Get q4
        help_q4 = atan2(-(w[3]*cos(s_q1) + w[4]*sin(s_q1)), -w[5])

        s_q4_1 = -s_q2_1 - s_q3_1 + help_q4
        s_q4_2 = -s_q2_2 - s_q3_2 + help_q4

        s_all = [
            np.array([self.wrapTau(s_q1), self.wrapTau(s_q2_1), self.wrapTau(s_q3_1), self.wrapTau(s_q4_1)]),
            np.array([self.wrapTau(s_q1), self.wrapTau(s_q2_2), self.wrapTau(s_q3_2), self.wrapTau(s_q4_2)]),
        ]
        
        return s_all

    def get_closest(self, q_all, q0):
        # Find closest IK solution to robot pose q0
        # INPUT (1): all IK sollutions, 6xN
        # INPUT (2): Current joint state configuration
        # OUTPUT: q 6x1
        q0 = np.array(q0);
        return min(q_all, key = lambda q: np.linalg.norm(q - q0));
    pass

    def wrapTau(self, x):
        return (x-tau*floor(x/(tau)+0.5))
    

    # LAB 2
    def taylor_path(self, w1, w2, q0, tol=0.01): 
        w1 = np.array(w1);
        w2 = np.array(w2);
        q0 = np.array(q0);
        # print(r(q0));
        # print(r(np.array([self.get_dk(q0), w1, w2])));
        # print();
        
        q1 = self.get_closest(self.get_ik(w1), q0);
        q2 = self.get_closest(self.get_ik(w2), q0);
        # this causes recursion overflow
        # return [*self.taylor_path_rec(w1, w2, q0, tol), q2];
        last_node = LinkedNode(q2);
        todo = [(q1, q2, w1, w2, last_node)];
        _count = 1;
        while todo:
            (q1, q2, w1, w2, node_after) = todo.pop();
            qm = (q1 + q2) / 2;
            wm_actual1= self.get_dk(q1);
            wm_actual = self.get_dk(qm);
            wm_actual2= self.get_dk(q2);
            wm_wanted = (w1 + w2) / 2;
            dist = np.linalg.norm(wm_actual - wm_wanted);
            if dist <= tol:
                continue;
            else:
                # _count += 1;
                # if _count > 10:
                #     print(dist);
                #     raise BaseException;
                # pass
                qm = self.get_closest(self.get_ik(wm_wanted), q1);
                # print(r(np.array([q1, q2, (q1 + q2) / 2, qm])), r(np.array([w1, w2, wm_actual1, wm_actual, wm_actual2, wm_wanted])), dist, sep = "\n", end = "\n" + "-"*50+"\n");
                node = node_after.insertBefore(qm);
                todo.append((q1, qm, w1, wm_wanted, node));
                todo.append((qm, q2, wm_wanted, w2, node_after));
            pass
        pass
        qs = [last_node.data];
        while last_node.prev:
            last_node = last_node.prev;
            qs.append(last_node.data);
        pass
        qs.reverse();
        return qs;
    pass
    def taylor_path_rec(self, w1, w2, q0, tol=0.01):
        # Funkcija implementira Taylorov postupak
        # ULAZI Funkcije:
        #   - w1: Kartezijska poza pocetne tocke, izrazena kao numpy vektor 6x1
        #   - w2: Kartezijska poza krajnje tocke, izrazena kao numpy vektor 6x1
        #   - q0: Pocetna poza zglobova robota, izrazena kao numpy vektor Nx1 
        #   - tol: Zadana tolerancija, izrazena kao Float64
        # IZLAZI Funkcije: Tocke putanje u prostoru zglobova, izrazene kao numpy matrica, 
        #                  gdje je svaki novi red nova tocka.

        # Odredi inverz
        w1 = np.array(w1);
        w2 = np.array(w2);
        q0 = np.array(q0);
        
        q1 = self.get_closest(self.get_ik(w1), q0);
        q2 = self.get_closest(self.get_ik(w2), q0);
        
        qm = (q1 + q2) / 2;
        wm_actual = self.get_dk(qm);
        wm_wanted = (w1 + w2) / 2;
        if (np.linalg.norm(wm_actual - wm_wanted) <= tol):
            return [];
        else:
            qm = self.get_closest(self.get_ik(wm_wanted), q1);
            return [
                *self.taylor_path_rec(w1, wm_wanted, q1, tol), 
                qm, 
                *self.taylor_path_rec(wm_wanted, w2, qm, tol),
            ];
        pass
    pass

    def interpolate_q(self, Q, T_param):
        # Polinomska interpolacija jedinstvenim polinomom u prostoru zglobova.
        # Svaki od 4 reda ulazne matrice Q predstavlja vektor vrijednosti zglobova
        # kroz koji manipulator mora proci. Izlaz funkcije su vrijednosti 
        # otipkanog polinoma frekvencijom 10 Hz.
        # ULAZI Funkcije:
        #   - Q: Točke putanje u prostoru zglobova, izrazena kao numpy matrica Nx6
        #   - T_param: Parametričko vrijeme segmenta
        # IZLAZI Funkcije: Matrica točaka zglobova, otipkanog polinoma 
        #       frekvencijom 10 Hz.
        global polynomes;
        polynomes = [];
        v = 0; a = 0;
        for (q1, q2, T) in zip(Q, Q[1 : ], T_param):
            P = Polynome([q1, v, a/2, (q2-q1-v*T-a*T**2/2)/T**3]);
            D = P.derivative();
            DD = D.derivative();
            polynomes.append(P);
            v = D(T);
            a = DD(T);
        pass
        q1 = q2;
        q2 = Q[0];
        T = T_param[-1];
        polynomes.append(Polynome([q1, v, a/2, 
            (+ 20*q2 - 20*q1 - 12*v*T - 3*a*T**2)/(2*T**3),
            (- 30*q2 + 30*q1 + 16*v*T + 3*a*T**2)/(2*T**4),
            (+ 12*q2 - 12*q1 -  6*v*T -   a*T**2)/(2*T**5),
        ]));
        
        if TESTING:
            assert close(polynomes[0 ].derivative()(0          ), np.zeros(4));
            assert close(polynomes[-1].derivative()(T_param[-1]), np.zeros(4));
            assert close(polynomes[0 ].derivative().derivative()(0          ), np.zeros(4));
            assert close(polynomes[-1].derivative().derivative()(T_param[-1]), np.zeros(4));
            for (q1, q2, P, T) in zip(Q, Q[1 : ] + [Q[0]], polynomes, T_param):
                assert close(P(0), q1), T;
                assert close(P(T), q2), T;
            pass
            for (P1, P2, T1) in zip(polynomes, polynomes[1 : ], T_param):
                assert close(P1.derivative()             (T1), P2.derivative()             (0)), T1;
                assert close(P1.derivative().derivative()(T1), P2.derivative().derivative()(0)), T1;
            pass
        pass
    
        q_all = [];
        t = 0;
        for (P, T) in zip(polynomes, T_param, strict = True):
            while t < T:
                q = P(t);
                q_all.append(q);
                t += 0.1;
            pass
            t -= T;
        pass
        return q_all;
    pass

    def maxLinear(self, P, T): return max(abs(P(0)), abs(P(T)));
    def maxQuadratic(self, P, T):
        ys = [abs(P(0)), abs(P(T))];
        xM = - P[1] / (2*P[2]);
        if 0 <= xM <= T: ys.append(abs(P(xM)));
        return max(ys);
    pass
    def maxCubic(self, P, T):
        ys = [abs(P(0)), abs(P(T))];
        x_mid = - P[2];
        x_off = sqrt(P[2]**2 - 3*P[3]*P[1]);
        xM1 = (x_mid - x_off) / (3*P[3]);
        xM2 = (x_mid + x_off) / (3*P[3]);
        if 0 <= xM1 <= T: ys.append(abs(P(xM1)));
        if 0 <= xM2 <= T: ys.append(abs(P(xM2)));
        return max(ys);
    pass

    def ho_cook(self, Q, v_max_lim, a_max_lim, f_s=250):
        # Funkcija implementira Ho-Cookovu metodu planiranja trajektorije
        # robotske ruke postivajuci zadana ogranicenja brzine (v_max) i akceleracije
        # (a_max). Funkcija kao rezultat vraca otipkanu trajektoriju frekvencijom
        # f_s u prostoru zglobova
        # ULAZI Funkcije:
        #   - Q: Točke putanje u prostoru zglobova, izrazena kao numpy matrica Nx6
        #   - v_max: Ograničenja brzina zglobova, izrazeno kao vektor 6x1
        #   - a_max: Ograničenja akceleracija zglobova, izrazeno kao vektor 6x1
        #   - f_s: Frekvencija otipkavanja trajektorije
        # IZLAZI Funkcije: Otipkane točke trajektorije, izrazene kao numpy matrica, 
        #                  gdje je svaki novi red nova tocka.
        
        Q = np.array(Q);
        v_max_lim= np.array(v_max_lim);
        a_max_lim= np.array(a_max_lim);

        (n, nq) = np.shape(Q);
        if TESTING: tprint(n, nq);

        # Racunanje parametrickog vremena (5.24)
        Ts = np.array([sqrt(sum((Q[j][i] - Q[j - 1][i])**2 for i in range(nq))) for j in range(1, n)]);

        # Prva iteracija
        iter_max = 10;
        iter_cnt = 0;
        S = 0

        # Postavi uvjet za loop petlju dok Ho-Cook nije zadovoljen
        for iter_cnt in range(iter_max):
            if TESTING: tprint(iter_cnt, "-" * 40);
            if TESTING: tprint(r(Ts));
            # Izracunaj matrie M, A i Dq (5.50)
            # Popuni matricu M
            MT = np.zeros((n - 2, n - 2));
            MT[0  ][0  ] = 3/Ts[0  ] + 2/Ts[1  ];
            MT[0  ][1  ] =             1/Ts[1  ];
            MT[n-3][n-3] = 3/Ts[n-2] + 2/Ts[n-3];
            MT[n-3][n-4] =             1/Ts[n-3];
            for i in range(1, n - 3):
                MT[i][i-1] = Ts[i + 1];
                MT[i][i  ] = 2 * (Ts[i + 1] + Ts[i]);
                MT[i][i+1] = Ts[i];
            pass
            if TESTING: tprint(r(np.transpose(MT)));

            # Popuni matricu A
            AT = np.array([
                + 3 / Ts[1]**2 * (Q[2] - Q[1])
                + 6 / Ts[0]**2 * (Q[1] - Q[0])
            ] + [
                3 / (Ts[i-1]*Ts[i]) * (
                    + Ts[i-1]**2 * (Q[i+1] - Q[i])
                    + Ts[i  ]**2 * (Q[i] - Q[i-1])
                )
                for i in range(2, n - 2)
            ] + [
                + 6 / Ts[n-2]**2 * (Q[n-1] - Q[n-2])
                + 3 / Ts[n-3]**2 * (Q[n-2] - Q[n-3])
            ]);
            if TESTING: tprint(r(np.transpose(AT)));

            # Izracunaj Dq
            DqT = np.linalg.inv(MT) @ AT;
            if TESTING: tprint(r(np.transpose(DqT)));


            # Izracunaj matricu B
            B = [None] * (n - 1); # could have [] with append, but whatever

            # Prvi segment (5.35)
            B[0] = np.transpose([Q[0], Q[1], DqT[0]]) @ np.array([
                [1, 0, 0, -4/Ts[0]**3, +3/Ts[0]**4],
                [0, 0, 0, +4/Ts[0]**3, -3/Ts[0]**4],
                [0, 0, 0, -1/Ts[0]**2, +1/Ts[0]**3],
            ]);
            if TESTING: tprint("B[0]");
            if TESTING: tprint(r(B[0]));

            # Medjusegmenti (5.23)
            for i in range(1, n - 2):
                B[i] = np.transpose([Q[i], Q[i+1], DqT[i-1], DqT[i]]) @ np.array([
                    [1, 0, -3/Ts[i]**2, +2/Ts[i]**3],
                    [0, 0, +3/Ts[i]**2, -2/Ts[i]**3],
                    [0, 1, -2/Ts[i]   , +1/Ts[i]**2],
                    [0, 0, -1/Ts[i]   , +1/Ts[i]**2],
                ]);
                if TESTING: tprint(f"B[{i}]");
                if TESTING: tprint(r(B[i]));
            pass


            # Zadnji segment (5.41)
            B[n - 2] = np.transpose([Q[n-2], Q[n-1], DqT[n - 3]]) @ np.array([
                [1, 0, -6/Ts[n-2]**2, +8/Ts[n-2]**3, -3/Ts[n-2]**4],
                [0, 0, +6/Ts[n-2]**2, -8/Ts[n-2]**3, +3/Ts[n-2]**4],
                [0, 1, -3/Ts[n-2]   , +3/Ts[n-2]**2, -1/Ts[n-2]**3],
            ]);
            if TESTING: tprint(f"B[{n - 2}]");
            if TESTING: tprint(r(B[n - 2]));
            B = [[Polynome(row) for row in b] for b in B];

            # Odredi max. brzine i akceleracije
            max_speeds = [];
            max_accels = [];
            Ds = [row.derivative() for row in B[0]];
            DDs = [D.derivative() for D in Ds];
            max_speeds.append([self.maxCubic    (D , Ts[0]) for D  in Ds ]);
            max_accels.append([self.maxQuadratic(DD, Ts[0]) for DD in DDs]);
            
            for i in range(1, n - 2):
                Ds = [row.derivative() for row in B[i]];
                DDs = [D.derivative() for D in Ds];
                max_speeds.append([self.maxQuadratic(D , Ts[i]) for D  in Ds ]);
                max_accels.append([self.maxLinear   (DD, Ts[i]) for DD in DDs]);
            pass
            
            Ds = [row.derivative() for row in B[n-2]];
            DDs = [D.derivative() for D in Ds];
            max_speeds.append([self.maxCubic    (D , Ts[n-2]) for D  in Ds ]);
            max_accels.append([self.maxQuadratic(DD, Ts[n-2]) for DD in DDs]);
            
            if TESTING: tprint();
            if TESTING: tprint(r(np.transpose(max_speeds)));
            if TESTING: tprint(r(np.transpose(max_accels)));
            
            max_speed = [max(mxv[i] for mxv in max_speeds) for i in range(nq)];
            max_accel = [max(mxa[i] for mxa in max_accels) for i in range(nq)];
            
            if TESTING: tprint();
            if TESTING: tprint(r(max_speed));
            if TESTING: tprint(r(max_accel));
            
            # Odredi parametre uvjeta Sv i Sa
            Sv =      max(max_speed / v_max_lim) ;
            Sa = sqrt(max(max_accel / a_max_lim));
            
            S = max(Sv, Sa);
            if TESTING: tprint(r(Sv), r(Sa), S, S - 1);
            if abs(S - 1) < 1e-14: break;
            Ts *= S;
        pass


        # Otipkavanje trajektorije
        dt = 1/f_s

        q_all = [];
        t = 0;
        for (b, T) in zip(B, Ts, strict = True):
            while t < T:
                q = [p(t) for p in b];
                q_all.append(q);
                t += dt;
            pass
            t -= T;
        pass
        return q_all;

        # Prvi segment

        # Medjusegmenti segment

        # Zadnji segment
        
        # Spoji segmente i vrati rezultat u obliku
        #return [Q_q, Q_dq, Q_ddq]
        pass 
    pass
pass
class ORLAB_OpenManipulator(RoboMath):
    def __init__(self):
        super().__init__();
        
        # Define publishers
        self.file_PTP = '/home/kukapc/kuka_ws/src/OR_lab2_kuka/roblab_kuka/scripts/lab{NAD}_ptp.txt'

        # Define subscribers

        # Define services  --- Comment this if working offline, without robot
        self.open_manipulator_send_command_service = rospy.ServiceProxy('/open_manipulator/goal_joint_space_path', SetJointPosition)
    pass
    # GENERAL FUNCTIONS
    def moveRobot(self, q, t):
        serviceReq = SetJointPositionRequest()

        print (serviceReq)

        serviceReq.joint_position.joint_name.append('joint1')
        serviceReq.joint_position.joint_name.append('joint2')
        serviceReq.joint_position.joint_name.append('joint3')
        serviceReq.joint_position.joint_name.append('joint4')
        serviceReq.joint_position.position = [q[0], q[1], q[2], q[3]]

        serviceReq.path_time = t

        self.open_manipulator_send_command_service.call(serviceReq)
    pass
pass

W = [
    ( 0.035, -0.035,  0.030, 0, 0, -1),
    ( 0.070, -0.035,  0.030, 0, 0, -1),
    ( 0.070,  0.035,  0.030, 0, 0, -1),
    ( 0.035,  0.035,  0.030, 0, 0, -1),
];
v_max = [0.7] * 4;
a_max = [1] * 4;

r = lambda a: np.round(a, 4)
if TESTING:
    import random;
    import traceback;
    x = RoboMath();
    
    def close(a, b):
        try: return all(close(aa, bb) for (aa, bb) in zip(a, b));
        except: return a - b < 1e-3;
    pass
    def distToLine(w1, w2, w):
        w1 = np.array(w1);
        w2 = np.array(w2);
        w  = np.array(w );
        v = w2 - w1;
        u = w  - w1;
        return np.linalg.norm(u - (np.linalg.vecdot(u, v) * (v / np.linalg.norm(v)**2)));
    pass
    for I in range(100):
        try:
            q0 = [x.wrapTau(random.uniform(0, tau)) for i in range(4)];
            w0 = x.get_dk(q0);
            qs = x.get_ik(w0);
            q0 = x.get_closest(qs, q0);
            w0 = x.get_dk(q0);
            qs = x.get_ik(w0);
            q1 = x.get_closest(qs, q0);
            assert close(q0, q1), (q0, q1);
            ws = [x.get_dk(q) for q in qs];
            for (i, w) in enumerate(ws): assert close(w, w0), (i, w, w0);
        except ValueError as err:
            if err.args == ("math domain error",): continue;
            raise;
        pass
    pass

    q_prev = [0, 0, 0, 0];
    Q = [];
    for w in W:
        q = x.get_closest(x.get_ik(w), q_prev);
        Q.append(q);
        q_prev = q;
    pass
    T_param = [round(max(abs(q2 - q1))/0.1) for (q1, q2) in zip(Q, Q[1 : ] + [Q[0]])];
    print(r(T_param));

    for (w1, w2, q1) in zip(W, W[1 : ] + W[0], Q):
        tol = 0.001;
        qt = x.taylor_path(w1, w2, q1, tol);
        for (i, q) in enumerate(qt): 
            w = x.get_dk(q);
            assert distToLine(w1, w2, w) <= tol, (i, q);
            if i == 0: continue;
            qq = qt[i - 1];
            ww = x.get_dk(qq);
            qm = (q + qq) / 2;
            wm = (w + ww) / 2;
            wq = x.get_dk(qm);
            assert np.linalg.norm(wm - wq) <= tol, (i, wm, wq);
            assert distToLine(w1, w2, wq) <= tol, (i, wq);
        pass
        print(r(qt));
    pass
    
    qp = x.interpolate_q(Q, T_param);
    
    class HoCook: # Ho-Cook
        print("=== Ho-Cook ====================================");
        Q = [
            ( 0     , 0     , 0      ),
            (-0.2414, 0.8957,   pi/12),
            (-0.23  , 1.231 ,   pi/ 6),
            (-0.1433, 1.4455,   pi/ 4),
            ( 0     , pi/2  ,   pi/ 3),
            ( 0.1847, 1.6125, 5*pi/12),
            ( 0.3948, pi/2  ,   pi/ 2),
        ];
        v_max = [1.7453, 1.5708, 1.3963];
        a_max = [0.7854, 0.8727, 0.6981];
        ch = x.ho_cook(Q, v_max, a_max);
    pass
pass


if __name__ == '__main__' and not TESTING:

    rospy.init_node('ROBLAB_open_manipulator')
    node = ORLAB_OpenManipulator()

    q_prev = [0, 0, 0, 0];
    Q = [];
    for w in W:
        q = node.get_closest(node.get_ik(w), q_prev);
        Q.append(q);
        q_prev = q;
    pass
    T_param = [round(max(abs(q2 - q1))/0.1) for (q1, q2) in zip(Q, Q[1 : ] + [Q[0]])];

    qt = [Q[0]];
    for (w1, w2, q1) in zip(W, W[1 : ] + W[0], Q):
        tol = 0.001;
        qs = node.taylor_path(w1, w2, q1, tol);
        qt.extend(qs);
    pass
    
    qp = node.interpolate_q(Q, T_param);
    
    qhc = node.ho_cook(qt, v_max, a_max, 250);
    
    for q in qt:
        node.moveRobot(q, 1); # one step every second
    pass
    for q in qp:
        node.moveRobot(q, 0.1); # 10 steps every second (the same period as used in interpolate_q)
    pass
    for q in qhc:
        node.moveRobot(q, 1/250); # 250 steps every second (same as in ho_cook)
    pass
pass
classdef MpcControl_x < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(H/Ts); % Horizon steps
            N = N_segs + 1;      % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            Q = eye(nx);
            R = 10 * eye(nu);

            % constraints of the input 
            M = [1;-1];
            m = [deg2rad(15); deg2rad(15)];

            % constraints of the states
            F = [0 1 0 0; 0 -1 0 0];
            f = [deg2rad(10); deg2rad(10)];

            [K,Qf,~] = dlqr(mpc.A, mpc.B, Q, R);
            K = -K;

            Xf = polytope([F; M*K], [f;m]);
            Acl = [mpc.A + mpc.B*K];
            while 1
                prevXf = Xf;
                [T, t] = double(Xf);
                preXf = polytope(T*Acl, t);
                Xf = intersect(Xf, preXf);
                if isequal(prevXf, Xf)
                    break
                end
            end

            [Ff, ff] = double(Xf);

            % Plot the terminal invarient set
            
            figure;
            
            subplot(2,3,1);
            Xf.projection([1,2]).plot();
            xlabel('w_y [rad/s]');
            ylabel('beta [rad]');
            
            subplot(2,3,2);
            Xf.projection([1, 3]).plot();
            xlabel('w_y [rad/s]');
            ylabel('v_x [m/s]');
            
            subplot(2,3,3);
            Xf.projection([1,4]).plot();
            xlabel('w_y [rad/s]');
            ylabel('x [m]');
            
            subplot(2,3,4);
            Xf.projection([2,3]).plot();
            xlabel('beta [rad]');
            ylabel('v_x [m/s]');
            
            subplot(2,3,5);
            Xf.projection([2,4]).plot();
            xlabel('beta [rad]');
            ylabel('x [m]');
            
            subplot(2,3,6);
            Xf.projection([3,4]).plot();
            xlabel('v_x [m/s]');
            ylabel('x [m]');
            
            sgtitle("Terminal invariant set of x");

            A = mpc.A;
            B = mpc.B;
        
            obj = U(:,1)'*R*U(:,1);
            con = (X(:,2) == A*X(:,1) + B*U(:,1)) + (M*U(:,1) <= m);
            for i = 2:N-1
                con = con + (X(:,i+1) == A*X(:,i) + B*U(:,i));
                con = con + (F*X(:,i) <= f) + (M*U(:,i) <= m);
                obj = obj + X(:,i)'*Q*X(:,i) + U(:,i)'*R*U(:,i);
            end 
            obj = obj + X(:,N)'*Qf*X(:,N);
            con = con + (Ff*X(:,N) <= ff);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref}, {U(:,1), X, U});
        end
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);

            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            obj = 0;
            con = [xs == 0, us == 0];
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end

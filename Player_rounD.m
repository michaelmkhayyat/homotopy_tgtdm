classdef Player_rounD < handle

    properties
        ID

        opt

        initialStateGlobal
        currentStateGlobal
        goalStateGlobal

        initialStateFrenet
        currentStateFrenet
        goalStateFrenet

        stateHistory
        controlHistory
        
        params              % These are the parameters of the vehicle in terms of dynamics, kinematics, as well as state limits.
        wayPoints           % These are the indicatory waypoints which define the path that the vehicle needs to take as it goes from initialPosition to goalPosition
        referencePath       % This is the reference path
        scenarioActor

        pathInfo


    end

    methods

        function obj = Player_rounD(ID, wayPoints,initialStateFrenet, goalStateFrenet, wid, len, hei, color,N,dT, meterPerPixel)

            obj.ID = ID;
            obj.opt.params.N = N;
            obj.opt.params.dT = dT;
            
            obj.params = [];
            obj.params.draw.meterPerPixel = meterPerPixel;

            obj.params.width = wid;
            obj.params.length = len;
            obj.params.height = hei;
            safetyFactor = 1.2;
            obj.params.safetyWidth = (safetyFactor-1)*len +  wid;
            obj.params.safetyLength = safetyFactor*len;
            obj.params.safetyCircleDiameter = sqrt(obj.params.safetyWidth^2 + obj.params.safetyLength^2);
            obj.params.col = color;

            obj.wayPoints = wayPoints;
            obj.initialStateFrenet = initialStateFrenet; %Given as s, s dot.
            obj.currentStateFrenet = initialStateFrenet;
            obj.goalStateFrenet = goalStateFrenet;
            
            % Generated a reference path by fitting a piecewise-continuous clothoid spline to the waypoints
            obj.referencePath = referencePathFrenet(obj.wayPoints);
            
            obj.initializeOptimizationAndParameterVariables()
            obj.initializeCostFunctionAndConstraints()

            obj.initialStateGlobal = frenet2global(obj.referencePath,[initialStateFrenet(1) initialStateFrenet(2) 0 0 0 0]);
            obj.currentStateGlobal = frenet2global(obj.referencePath,[obj.currentStateFrenet(1) obj.currentStateFrenet(2) 0 0 0 0]);
            obj.goalStateGlobal = frenet2global(obj.referencePath,[goalStateFrenet(1) goalStateFrenet(2) 0 0 0 0]);
            
            lengthPath = obj.referencePath.PathLength;
            discretizationDistance = obj.referencePath.DiscretizationDistance;
            interpolationPoints = 0:discretizationDistance:lengthPath;

            interpolatedWayPoints = [];
            for i = 1:1:length(interpolationPoints)
                temp = interpolate(obj.referencePath, interpolationPoints(i));
                interpolatedWayPoints = [interpolatedWayPoints; temp(1:2)];
            end
            
            obj.pathInfo.centerPoints = interpolatedWayPoints;
            calculateExactPathBound = 1;

            if calculateExactPathBound == 1

                obj.calculateUpperAndLowerBound()

            else

                s1_upper = [];
                s1_lower = [];
    
                for i = 1:1:length(interpolationPoints)
                    temp1_l = frenet2global(obj.referencePath,[interpolationPoints(i) 0 0 -obj.params.safetyWidth/2 0 0]);
                    temp1_h = frenet2global(obj.referencePath,[interpolationPoints(i) 0 0 obj.params.safetyWidth/2 0 0]);
    
                    s1_upper = [s1_upper; temp1_h];
                    s1_lower = [s1_lower; temp1_l];
    
                end
    
                 obj.pathInfo.upperBound = s1_upper;
                 obj.pathInfo.lowerBound = s1_lower;

            end

        end

        function [] = initializeOptimizationAndParameterVariables(Player_rounD)

            import casadi.*

            n_states = 2;
            n_control = 1;

            x = sdpvar(Player_rounD.opt.params.N+1,n_states);    % States vector for Player_rounD 
            u = sdpvar(Player_rounD.opt.params.N,n_control);     % Control vector for Player_rounD 
            p = sdpvar(1,n_states);                        % Parameter vector for vector for Player_rounD 


            Player_rounD.opt.vars.x = x;
            Player_rounD.opt.vars.u = u;
            Player_rounD.opt.vars.p = p;

            Player_rounD.opt.discreteflag.x = zeros(Player_rounD.opt.params.N+1,n_states);
            Player_rounD.opt.discreteflag.u = zeros(Player_rounD.opt.params.N,n_control);
            Player_rounD.opt.discreteflag.p = zeros(1, n_states);

        end

        function [] = initializeCostFunctionAndConstraints(Player_rounD)
            
            import casadi.*

            objective = 0;
            P = 1;
            r = 5;

            vmax = 15; % Maximum Velocity
            amax  = 10; %Maximum Acceleration
            smax = Player_rounD.referencePath.PathLength;

            x = Player_rounD.opt.vars.x;
            u = Player_rounD.opt.vars.u;
            p = Player_rounD.opt.vars.p;

            for k = 1:1:Player_rounD.opt.params.N
                objective = objective + P*u(k,1)*u(k,1);
            end
            objective = objective - r*(x(k+1,1) - x(1,1));

            upperStateBound = repmat([smax vmax], Player_rounD.opt.params.N+1, 1);
            lowerStateBound = repmat([0 0], Player_rounD.opt.params.N+1, 1);
            upperControlBound = repmat([amax], Player_rounD.opt.params.N, 1);
            lowerControlBound = repmat([-amax], Player_rounD.opt.params.N, 1);

            Player_rounD.opt.costFunction = objective;
            Player_rounD.opt.constraints.state.vars = x(:);
            Player_rounD.opt.constraints.state.ub = upperStateBound(:);
            Player_rounD.opt.constraints.state.lb = lowerStateBound(:);    
            Player_rounD.opt.constraints.control.vars = u(:);
            Player_rounD.opt.constraints.control.ub = upperControlBound(:);
            Player_rounD.opt.constraints.control.lb = lowerControlBound(:);

            
            temp1 = [];
            temp2 = [];
            temp3 = [];
            for k = 1:1:Player_rounD.opt.params.N+1

                temp1 = [temp1; x(k,2)^2 + 2*amax*x(k,1)];
                temp2 = [temp2; 0];
                temp3 = [temp3; 2*amax*smax];
            
            end
            
            Player_rounD.opt.recursiveFeasabilityConstraints.vars = temp1;
            Player_rounD.opt.recursiveFeasabilityConstraints.lb = temp2;
            Player_rounD.opt.recursiveFeasabilityConstraints.ub = temp3;

            dynamics = [];

            useForwardEuler = 1;

            A_cont = [0 1; 0 0];
            b_cont = [0;1];

            if useForwardEuler == 1

                for k = 1:1:Player_rounD.opt.params.N
                   dynamics = [dynamics; x(k+1,1) - x(k,1) - Player_rounD.opt.params.dT*x(k,2) - Player_rounD.opt.params.dT^2/2*u(k,1)];
                   dynamics = [dynamics; x(k+1,2) - x(k,2) - Player_rounD.opt.params.dT*u(k,1)];
                end
                
                dynamics = [dynamics; x(1,1) - p(1)];
                dynamics = [dynamics; x(1,2) - p(2)];

            else

                for k = 1:1:Player_rounD.opt.params.N
                    k1 = A_cont * (x(k,1:2)') + b_cont * u(k,1);
                    k2 = A_cont * (x(k,1:2)' + Player_rounD.opt.params.dT/2 * k1) + b_cont * u(k,1);
                    k3 = A_cont * (x(k,1:2)' + Player_rounD.opt.params.dT/2 * k2) + b_cont * u(k,1);
                    k4 = A_cont * (x(k,1:2)' + Player_rounD.opt.params.dT*k3) + b_cont * u(k,1);
                    dynamics = [dynamics; x(k+1,1:2)' - x(k,1:2)' - Player_rounD.opt.params.dT/6 *( k1 + 2*k2 + 2*k3 + k4)];
 
                end
                
                dynamics = [dynamics; x(1,1) - p(1)];
                dynamics = [dynamics; x(1,2) - p(2)];

            end

            Player_rounD.opt.constraints.dynamics = dynamics;

        end
        
        function [] = plotReferencePath(Player_rounD,parent)
            hold on
            if nargin == 2
                plot(Player_rounD.pathInfo.centerPoints(:,1), Player_rounD.pathInfo.centerPoints(:,2), "Color", Player_rounD.params.col, "LineWidth", 1.5, "Parent",parent)
            else
                plot(Player_rounD.pathInfo.centerPoints(:,1), Player_rounD.pathInfo.centerPoints(:,2), "Color", Player_rounD.params.col,"LineWidth", 1.5)
            end
        end

        function[] = plotBounds(Player_rounD, parent)
            hold on
            plot(Player_rounD.pathInfo.upperBound(:,1), Player_rounD.pathInfo.upperBound(:,2), 'Color', [Player_rounD.params.col,0.5], 'LineWidth', 1, 'Parent',parent);
            plot(Player_rounD.pathInfo.lowerBound(:,1), Player_rounD.pathInfo.lowerBound(:,2), 'Color', [Player_rounD.params.col,0.5], 'LineWidth', 1, 'Parent',parent);
        end

        function [] = calculateUpperAndLowerBound (Player_rounD)
            
            lengthPath = Player_rounD.referencePath.PathLength;
            discretizationDistance = Player_rounD.referencePath.DiscretizationDistance;
            interpolationPoints = 0:discretizationDistance:lengthPath;

            Edge1 = zeros(length(interpolationPoints),2);
            Edge2 = zeros(length(interpolationPoints),2);
            Edge3 = zeros(length(interpolationPoints),2);
            Edge4 = zeros(length(interpolationPoints),2);

            edge1 =[];
            edge2 =[];
            edge3 =[];
            edge4 =[];
        
            for kk=1:1:length(interpolationPoints)
                state = frenet2global(Player_rounD.referencePath,[interpolationPoints(kk) 0 0 0 0 0]);
                temp = Player_rounD.getVertices([state(1), state(2), state(3)]);
                edge1 = [edge1 ; temp(1,:)];
                edge2  = [edge2 ; temp(2,:)];
                edge3 = [edge3 ; temp(3,:)];
                edge4  = [edge4 ; temp(4,:)];
            end

            edgesAlongPath = [edge1;edge2;edge3;edge4];
            k = boundary(edgesAlongPath(:,1), edgesAlongPath(:,2),1);

            enclosureX = edgesAlongPath(k,1);
            enclosureY = edgesAlongPath(k,2);
            s_positiveIndices = [];
            s_negativeIndices = [];
            s_arrangementPositive = [];
            s_arrangementNegative = [];

            for i = 1:1:length(enclosureX)
                temp = global2frenet(Player_rounD.referencePath, [enclosureX(i) enclosureY(i) 0 0 0 0]);

                if temp(4) > 0
                    s_arrangementPositive = [s_arrangementPositive; temp(1)];
                    s_positiveIndices = [s_positiveIndices; i];
                else
                    s_arrangementNegative = [s_arrangementNegative; temp(1)];
                    s_negativeIndices = [s_negativeIndices; i];
                end

            end

            [s_arrangementPositive,sortIndex] = sort(s_arrangementPositive);

            s_positiveIndices = s_positiveIndices(sortIndex);

            [s_arrangementNegative,sortIndex] = sort(s_arrangementNegative);

            s_negativeIndices = s_negativeIndices(sortIndex);
            
            Player_rounD.pathInfo.upperBound = [enclosureX(s_positiveIndices), enclosureY(s_positiveIndices)];
            Player_rounD.pathInfo.lowerBound = [enclosureX(s_negativeIndices), enclosureY(s_negativeIndices)];

            % 
%             leftEdgesAlongPath = [edge1;edge2];
%             rightEdgesAlongPath = [edge3;edge4];
%             k_left = boundary(leftEdgesAlongPath(:,1), leftEdgesAlongPath(:,2),1);
%             k_right = boundary(rightEdgesAlongPath(:,1), rightEdgesAlongPath(:,2),1);
% 
%             Player_rounD.pathInfo.upperBound = [leftEdgesAlongPath(k_left); rightEdgesAlongPath(k_left)];
%             Player_rounD.pathInfo.lowerBound = [leftEdgesAlongPath(k_right); rightEdgesAlongPath(k_right)];

        end


        function [h] = drawPlayer(Player_rounD, states) 

        H = Player_rounD.params.width/Player_rounD.params.draw.meterPerPixel;
        L = Player_rounD.params.length/Player_rounD.params.draw.meterPerPixel;
        xc=states(1);
        yc=states(2);
        theta= states(3);
        transformationMatrix= ([cos(theta), -sin(theta); sin(theta), cos(theta)]);
        X=([-L/2, L/2, L/2, -L/2]);
        Y=([-H/2, -H/2, H/2, H/2]);
        T = zeros(2,4);
        for i=1:4
            T(:,i)=transformationMatrix*[X(i); Y(i)];
        end
        
        x_lower_left=xc+T(1,1);
        x_lower_right=xc+T(1,2);
        x_upper_right=xc+T(1,3);
        x_upper_left=xc+T(1,4);
        y_lower_left=yc+T(2,1);
        y_lower_right=yc+T(2,2);
        y_upper_right=yc+T(2,3);
        y_upper_left=yc+T(2,4);
        x_coor=[x_lower_left x_lower_right x_upper_right x_upper_left];
        y_coor=[y_lower_left y_lower_right y_upper_right y_upper_left];
        triangleXPosition = [(x_coor(1)+x_coor(2))/2 (x_coor(3)+x_coor(4))/2 (x_coor(2)+x_coor(3))/2];
        triangleYPosition = [(y_coor(1)+y_coor(2))/2 (y_coor(3)+y_coor(4))/2 (y_coor(2)+y_coor(3))/2];
        h = patch('Vertices',[x_coor; y_coor]','Faces',[1 2 3 4],'Edgecolor','none'	,'Facecolor',Player_rounD.params.col,'FaceAlpha',1,'PickableParts', 'visible');
        hold on
        triangle = patch(triangleXPosition, triangleYPosition, [0 0 0],'LineJoin','round');
        end
        
        function vertices = getVertices(Player_rounD, stateVector)

            cx = stateVector(1);
            cy = stateVector(2);
            theta = stateVector(3);

            w = Player_rounD.params.safetyWidth;
            h =  Player_rounD.params.safetyLength;
            
            vertices = [cx - (h*cos(theta))/2 - (w*sin(theta))/2,...
                cx - (h*cos(theta))/2 + (w*sin(theta))/2,...
                cx + (h*cos(theta))/2 + (w*sin(theta))/2,...
                cx + (h*cos(theta))/2 - (w*sin(theta))/2;
                cy + (w*cos(theta))/2 - (h*sin(theta))/2,...
                cy - (w*cos(theta))/2 - (h*sin(theta))/2,...
                cy - (w*cos(theta))/2 + (h*sin(theta))/2,...
                cy + (w*cos(theta))/2 + (h*sin(theta))/2];
            
            vertices = vertices';
        end
    end

end
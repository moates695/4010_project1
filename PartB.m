% Marcus Oates 
% z5257541
% Date: 22/04/23
% File contains Part B code
% usage: >>Part_B('DataUsr_p020')

function main(file)
    load(file); 
    extract(data);
end

% ----------------------------------------

function extract(data)
    disp('Begin optimizing');
    gyroBias = optimize(data);
    disp('End optimizing');

    figure(11); clf();
    
    landmarks = data.Context.Landmarks;
    plot(landmarks(1,:), landmarks(2,:), 'ko');
    title('Global CF (Displaying Data)');
    xlabel('x (m)'); 
    ylabel('y (m)');

    hold on;
    walls = data.Context.Walls;
    plot(walls(1,:), walls(2,:), 'color', [0,1,0]*0.7, 'linewidth', 3);    
    
    p0=data.pose0;
    plot(p0(1),p0(2),'r*','markersize',10);

    plot(nan, nan,'g.');
    plot(nan, nan, 'r.');
    plot(nan, nan, 'b.');

    legend({'landmarks','walls (middle planes)','initial position', 'prediction (bias=0)', 'prediction (bias)', 'ground truth'});
    hold off;

    X1 = data.pose0; % [meters; meters; radians] with estimated bias
    vw1 = [0; 0]; % [meters/sec; rad/sec] with estimated bias
    X2 = data.pose0; % [meters; meters; radians] assuming bias = 0
    vw2 = [0; 0]; % [meters/sec; rad/sec] assuming bias = 0
    events  = data.table;
    event0 = events(:,1);
    t_last = 0.0001 * double(event0(1));    
    ground = data.verify.poseL;
        
    disp('Begin sampling events'); 
    subsampleIdx = 1;

    for i = 1:data.n
        event = events(:,i);
        index = event(2);
        sensorID = event(3);

        t_curr = 0.0001 * double(event(1));
        dt = t_curr - t_last;
        t_last = 0.0001 * double(event(1));

        X1 = kinematicModel(X1, vw1, dt);
        X2 = kinematicModel(X2, vw2, dt);

        switch sensorID
            case 1
                figure(11); hold on; legend('AutoUpdate','off');
                plot(X2(1), X2(2), 'g.');
                plot(ground(1,subsampleIdx), ground(2,subsampleIdx), 'b.');
                plot(X1(1), X1(2), 'r.');
                hold off;
                subsampleIdx = subsampleIdx + 1;
            case 2
                vw1 = data.vw(:, index);
                vw1(2) = vw1(2) + gyroBias;
                vw2 = data.vw(:, index);
            otherwise
                continue
        end
    end

    disp('End sampling events');
end

% --------------------------------------------------------------------------------

function X = kinematicModel(X, vw, dt)
    dXdt = [vw(1) * cos(X(3)); vw(1) * sin(X(3)); vw(2)];
    X = X + dXdt * dt;
end

function bias = optimize(data)
    thresh = 0.1*pi/180;
    low = -2*pi/180;
    high = 2*pi/180;
    bias = 0;
    ground = data.verify.poseL;
    events  = data.table;
    buffNum = 0;
    for i = 1:data.n
        event = events(:,i);
        sensorID = event(3);
        if sensorID == 1
            buffNum = buffNum + 1;
        end
    end

    lastBias = 100;

    while 1
        if abs(lastBias - bias) < 0.00001
            break
        end
        fprintf('Current bias estimation: %.8f deg/sec\n', -1*bias*180/pi);
        refresh = false;

        X = data.pose0;
        vw = [0; 0];

        events  = data.table;
        event0 = events(:,1);
        t_last = 0.0001 * double(event0(1)); 
        subsampleIdx = 1;
        %headings = zeros(buffNum, 1);
        headings = [];

        for i = 1:data.n
            event = events(:,i);
            index = event(2);
            sensorID = event(3);

            t_curr = 0.0001 * double(event(1));
            dt = t_curr - t_last;
            t_last = 0.0001 * double(event(1));

            X = kinematicModel(X, vw, dt);        
            switch sensorID 
                case 1
                    headings(end+1) = X(3);
                    if subsampleIdx < buffNum
                        subsampleIdx = subsampleIdx + 1;
                        continue
                    end
                    avg = sum(headings) / length(headings);
                    groundAvg = sum(ground(3,(subsampleIdx-buffNum+1):subsampleIdx)) / buffNum;
                    if angleBetween(avg, groundAvg) >= thresh
                        refresh = true;
                    else
                        subsampleIdx = subsampleIdx + 1;
                    end
                case 2
                    vw = data.vw(:, index);
                    vw(2) = vw(2) + bias;

                otherwise
                    continue
            end
            if refresh
                break
            end
        end
        if refresh
            if spinDirection(avg, groundAvg) == 1
                high = bias;
            else
                low = bias; 
            end
            lastBias = bias;
            bias = (low+high)/2;
            continue
        else
            break
        end
    end
    fprintf('The bias for the gyroscope is %.8f degrees/sec\n', -1*bias*180/pi);
end

function ang = angleBetween(a, b)
    ang = abs(a - b);
    if ang > pi/2
        ang = 2*pi - ang;
    end
end

function dir = spinDirection(a, b)
    if a < b
        if abs(a-b) <= 180
            dir = -1;
        else
            dir = 1;
        end
    else
        if abs(a-b) > 180
            dir = -1;
        else
            dir = 1;
        end
    end
end

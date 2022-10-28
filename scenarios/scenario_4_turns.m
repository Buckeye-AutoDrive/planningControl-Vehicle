function [allData, scenario, sensors] = scenario_4_turns()
%scenario_4_turns - Returns sensor detections
%    allData = scenario_4_turns returns sensor detections in a structure
%    with time for an internally defined scenario and sensor suite.
%
%    [allData, scenario, sensors] = scenario_4_turns optionally returns
%    the drivingScenario and detection generator objects.

% Generated by MATLAB(R) 9.11 (R2021b) and Automated Driving Toolbox 3.4 (R2021b).
% Generated on: 15-May-2022 01:54:49

% Create the drivingScenario object and ego car
[scenario, egoVehicle] = createDrivingScenario;

% Create all the sensors
[sensors, numSensors] = createSensors(scenario);

allData = struct('Time', {}, 'ActorPoses', {}, 'ObjectDetections', {}, 'LaneDetections', {}, 'PointClouds', {}, 'INSMeasurements', {});
running = true;
while running

    % Generate the target poses of all actors relative to the ego vehicle
    poses = targetPoses(egoVehicle);
    time  = scenario.SimulationTime;

    objectDetections = {};
    laneDetections   = [];
    ptClouds = {};
    insMeas = {};
    isValidTime = false(1, numSensors);

    % Generate detections for each sensor
    for sensorIndex = 1:numSensors
        sensor = sensors{sensorIndex};
        % Generate the ego vehicle lane boundaries
        if isa(sensor, 'visionDetectionGenerator')
            maxLaneDetectionRange = min(500,sensor.MaxRange);
            lanes = laneBoundaries(egoVehicle, 'XDistance', linspace(-maxLaneDetectionRange, maxLaneDetectionRange, 101));
        end
        [objectDets, numObjects, isValidTime(sensorIndex)] = sensor(poses, time);
        objectDetections = [objectDetections; objectDets(1:numObjects)]; %#ok<AGROW>
    end

    % Aggregate all detections into a structure for later use
    if any(isValidTime)
        allData(end + 1) = struct( ...
            'Time',       scenario.SimulationTime, ...
            'ActorPoses', actorPoses(scenario), ...
            'ObjectDetections', {objectDetections}, ...
            'LaneDetections', {laneDetections}, ...
            'PointClouds',   {ptClouds}, ... %#ok<AGROW>
            'INSMeasurements',   {insMeas}); %#ok<AGROW>
    end

    % Advance the scenario one time step and exit the loop if the scenario is complete
    running = advance(scenario);
end

% Restart the driving scenario to return the actors to their initial positions.
restart(scenario);

% Release all the sensor objects so they can be used again.
for sensorIndex = 1:numSensors
    release(sensors{sensorIndex});
end

%%%%%%%%%%%%%%%%%%%%
% Helper functions %
%%%%%%%%%%%%%%%%%%%%

% Units used in createSensors and createDrivingScenario
% Distance/Position - meters
% Speed             - meters/second
% Angles            - degrees
% RCS Pattern       - dBsm

function [sensors, numSensors] = createSensors(scenario)
% createSensors Returns all sensor objects to generate detections

% Assign into each sensor the physical and radar profiles for all actors
profiles = actorProfiles(scenario);
sensors{1} = visionDetectionGenerator('SensorIndex', 1, ...
    'MinObjectImageSize', [5 5], ...
    'DetectorOutput', 'Objects only', ...
    'Intrinsics', cameraIntrinsics([800 800],[320 240],[480 640]), ...
    'ActorProfiles', profiles);
sensors{2} = drivingRadarDataGenerator('SensorIndex', 2, ...
    'MountingLocation', [3.7 0 0.5], ...
    'RangeLimits', [0 100], ...
    'HasNoise', false, ...
    'TargetReportFormat', 'Detections', ...
    'HasElevation', true, ...
    'HasOcclusion', false, ...
    'HasFalseAlarms', false, ...
    'FieldOfView', [35 50], ...
    'Profiles', profiles);
numSensors = 2;

function [scenario, egoVehicle] = createDrivingScenario
% createDrivingScenario Returns the drivingScenario defined in the Designer

% Construct a drivingScenario object.
scenario = drivingScenario;

% Add all road segments
roadCenters = [-6.7 2.5 0;
    83 2.5 0];
marking = laneMarking('Unmarked');
laneSpecification = lanespec(1, 'Width', 6, 'Marking', marking);
road(scenario, roadCenters, 'Lanes', laneSpecification, 'Name', 'Road');

roadCenters = [23.4 3.2 0;
    23.5 -20.5 0];
marking = laneMarking('Unmarked');
laneSpecification = lanespec(1, 'Width', 6, 'Marking', marking);
road(scenario, roadCenters, 'Lanes', laneSpecification, 'Name', 'Road1');

roadCenters = [42.9 26.5 0;
    42.9 3 0];
marking = laneMarking('Unmarked');
laneSpecification = lanespec(1, 'Width', 6, 'Marking', marking);
road(scenario, roadCenters, 'Lanes', laneSpecification, 'Name', 'Road2');

roadCenters = [39.8 29.3 0;
    68.5 29.3 0];
marking = laneMarking('Unmarked');
laneSpecification = lanespec(1, 'Width', 6, 'Marking', marking);
road(scenario, roadCenters, 'Lanes', laneSpecification, 'Name', 'Road');

roadCenters = [68.7 32.4 0;
    68.7 -0.5 0];
marking = laneMarking('Unmarked');
laneSpecification = lanespec(1, 'Width', 6, 'Marking', marking);
road(scenario, roadCenters, 'Lanes', laneSpecification, 'Name', 'Road1');

% Add the ego vehicle
egoVehicle = vehicle(scenario, ...
    'ClassID', 1, ...
    'Position', [-5 1.4 0], ...
    'FrontOverhang', 0.9, ...
    'Wheelbase', 2.8, ...
    'Mesh', driving.scenario.carMesh, ...
    'PlotColor', [0 114 189] / 255, ...
    'Name', 'Car');
waypoints = [-5 1.4 0;
    9.4 1.4 0;
    19.9 1.4 0;
    28.9 1.4 0;
    36.5 1.4 0;
    40.5 1.4 0;
    43.39 2.08 0;
    44.34 6.1 0;
    44.2 11.1 0;
    44.2 17.1 0;
    44.2 23.1 0;
    45.33 27.62 0;
    48.6 28.2 0;
    52.6 28.2 0;
    58.6 28.2 0;
    62.6 28.2 0;
    66.53 27.39 0;
    67.7 23.5 0;
    67.7 18.5 0;
    67.7 14.5 0;
    67.7 10.5 0;
    67.7 6.5 0;
    68.66 2.59 0.01;
    72.6 1.5 0;
    76.6 1.5 0;
    80.6 1.5 0];
speed = [30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30;30];
waittime = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
trajectory(egoVehicle, waypoints, speed, waittime);

% Add the non-ego actors
actor(scenario, ...
    'ClassID', 7, ...
    'Length', 0.7, ...
    'Width', 0.2, ...
    'Height', 0.1, ...
    'Position', [26.4 1 5.1], ...
    'Pitch', 90, ...
    'PlotColor', [0 0 0] / 255, ...
    'Name', 'Traffic_Light');

actor(scenario, ...
    'ClassID', 7, ...
    'Length', 0.7, ...
    'Width', 0.2, ...
    'Height', 0.1, ...
    'Position', [81.9466666666667 1 5.1], ...
    'Pitch', 90, ...
    'PlotColor', [0 0 0] / 255, ...
    'Name', 'Traffic_Light');

actor(scenario, ...
    'ClassID', 7, ...
    'Length', 0.7, ...
    'Width', 0.2, ...
    'Height', 0.1, ...
    'Position', [69.9033333333334 27.85 5.1], ...
    'Pitch', 90, ...
    'PlotColor', [0 0 0] / 255, ...
    'Name', 'Traffic_Light');

actor(scenario, ...
    'ClassID', 7, ...
    'Length', 0.7, ...
    'Width', 0.2, ...
    'Height', 0.1, ...
    'Position', [45.28 1 5.1], ...
    'Pitch', 90, ...
    'PlotColor', [0 0 0] / 255, ...
    'Name', 'Traffic_Light');

actor(scenario, ...
    'ClassID', 3, ...
    'Length', 1.7, ...
    'Width', 0.45, ...
    'Height', 1.7, ...
    'Position', [10.3 4.1 0], ...
    'Mesh', driving.scenario.bicycleMesh, ...
    'PlotColor', [237 177 32] / 255, ...
    'Name', 'Bicycle');

actor(scenario, ...
    'ClassID', 8, ...
    'Length', 0.01, ...
    'Width', 0.7, ...
    'Height', 0.7, ...
    'Position', [26.79 -0.96 2.4], ...
    'PlotColor', [255 0 0] / 255, ...
    'Name', 'Stop Sign');

actor(scenario, ...
    'ClassID', 8, ...
    'Length', 0.01, ...
    'Width', 0.7, ...
    'Height', 0.7, ...
    'Position', [45.4 -1.02 2.4], ...
    'PlotColor', [255 0 0] / 255, ...
    'Name', 'Stop Sign');

actor(scenario, ...
    'ClassID', 9, ...
    'Length', 0.1, ...
    'Width', 0.1, ...
    'Height', 2.4, ...
    'Position', [26.79 -0.96 0], ...
    'PlotColor', [128 128 128] / 255, ...
    'Name', 'Pole');

actor(scenario, ...
    'ClassID', 9, ...
    'Length', 0.1, ...
    'Width', 0.1, ...
    'Height', 2.4, ...
    'Position', [45.395 -1.017 0], ...
    'PlotColor', [128 128 128] / 255, ...
    'Name', 'Pole1');

actor(scenario, ...
    'ClassID', 8, ...
    'Length', 0.01, ...
    'Width', 0.7, ...
    'Height', 0.7, ...
    'Position', [44.5 35 0], ...
    'Yaw', 90, ...
    'PlotColor', [255 0 0] / 255, ...
    'Name', 'Stop Sign');

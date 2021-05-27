%%
robot = importrobot('urdf/irb120_3_58.urdf');

[timeP1, forceData1, timeF1, positionData1]= read(1,robot);
[timeP2, forceData2, timeF2, positionData2]= read(2,robot);
[timeP3, forceData3, timeF3, positionData3]= read(3,robot);
[timeP4, forceData4, timeF4, positionData4]= read(4,robot);
[timeP5, forceData5, timeF5, positionData5]= read(5,robot);
%% Trial 1
zeroF1 = forceData1(1,3);
force1 = forceData1(:,3)-zeroF1;
force1 = -force1(1160:2321);
zeroP1 = positionData1(1160,3);
depth1 = positionData1(1160:2321,3)-zeroP1;
depth1 = -depth1;


%% Trial 2
zeroF2 = forceData2(1,3);
force2 = forceData2(:,3)-zeroF2;
force2 = -force2(2500:3800);
zeroP2 = positionData2(2500,3);
depth2 = positionData2(2500:3800,3)-zeroP2;
depth2 = -depth2;

%% Trial 3
zeroF3 = forceData3(1,3);
force3 = forceData3(:,3)-zeroF3;
force3 = -force3(2101:3301);
zeroP3 = positionData3(2101,3);
depth3 = positionData3(2101:3301,3)-zeroP3;
depth3 = -depth3;

%% Trial 4
zeroF4 = forceData4(1,3);
force4 = forceData4(:,3)-zeroF4;
force4 = -force4(2453:3730);
zeroP4 = positionData4(2453,3);
depth4 = positionData4(2453:3730,3)-zeroP4;
depth4 = -depth4;

%% Trial 5
zeroF5 = forceData5(1,3);
force5 = forceData5(:,3)-zeroF3;
force5 = -force5(2101:3361);
zeroP5 = positionData5(2101,3);
depth5 = positionData5(2101:3361,3)-zeroP5;
depth5 = -depth5;

figure

plot(depth1,force1)
hold on
plot(depth2,force2)
plot(depth3,force3)
plot(depth4,force4)
plot(depth5,force5)
legend('trial 1','trial 2','trial 3','trial 4','trial 5');

figure
plot(timeP2,positionData2(:,3))
legend('end effector z')

figure
plot(timeF2,forceData2)
legend('force x', 'force y', 'force z')


function [time1Data, forceData, time2Data, positionData]= read(trial,robot)

switch trial
    case 1
        %trial 1
        bag = rosbag('sanddata/2021-05-25-16-56-22.bag');
    case 2
        %trial 2
        bag = rosbag('sanddata/2021-05-25-17-03-53.bag');
    case 3
        %trial 3
        bag = rosbag('sanddata/2021-05-25-17-10-19.bag');
    case 4
        %trial 4
        bag = rosbag('sanddata/2021-05-25-17-24-57.bag');
    case 5
        %trial 5
        bag = rosbag('sanddata/2021-05-25-17-27-50.bag');    
end

startT = bag.StartTime;
endT = bag.EndTime;

%% Force sensor data

bagselectData = select(bag, 'Time', [startT endT], 'Topic', '/data');

ts2 = timeseries(bagselectData, 'Wrench.Force.X', 'Wrench.Force.Y','Wrench.Force.Z');
forceData = ts2.Data;
time2Data = ts2.Time-startT;

%% Joint states data

bagselectJointStates = select(bag, 'Time', [startT endT],'Topic', '/joint_states');
msgs1 = readMessages(bagselectJointStates);
numMess = bagselectJointStates.NumMessages;

positionData = zeros(numMess,3);
config = homeConfiguration(robot);

for i = 1:numMess
    msgstemp = msgs1{i};  
    config(1).JointPosition = msgstemp.Position(1);
    config(2).JointPosition = msgstemp.Position(2);
    config(3).JointPosition = msgstemp.Position(3);
    config(4).JointPosition = msgstemp.Position(4);
    config(5).JointPosition = msgstemp.Position(5);
    config(6).JointPosition = msgstemp.Position(6);
    T = getTransform(robot, config, 'tool0');
    positionData(i,:) = [T(1,end), T(2,end), T(3,end)];
end


ts1 = timeseries(bagselectJointStates);
time1Data = ts1.Time-startT;


end

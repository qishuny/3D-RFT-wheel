[time1Data, forceData, time2Data, positionData]= read(2);


figure
plot(time1Data,positionData(:,3))
legend('end effector z')


figure
plot(time2Data,forceData)
legend('force x', 'force y', 'force z')


function [time1Data, forceData, time2Data, positionData]= read(trial)

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

bagselect2 = select(bag, 'Time', [startT endT], 'Topic', '/data');
msgs2 = readMessages(bagselect2);
temp2 = msgs2{2};

ts2 = timeseries(bagselect2, 'Wrench.Force.X', 'Wrench.Force.Y','Wrench.Force.Z');
forceData = ts2.Data;
time2Data = ts2.Time-startT;

%%
robot = importrobot('urdf/irb120_3_58.urdf');

bagselect1 = select(bag, 'Time', [startT endT],'Topic', '/joint_states');
msgs1 = readMessages(bagselect1);
numMess = bagselect1.NumMessages;

positionData = zeros(numMess,3);
config = homeConfiguration(robot);

for i = 1:numMess
    temp1 = msgs1{i};
    
    config(1).JointPosition = temp1.Position(1);
    config(2).JointPosition = temp1.Position(2);
    config(3).JointPosition = temp1.Position(3);
    config(4).JointPosition = temp1.Position(4);
    config(5).JointPosition = temp1.Position(5);
    config(6).JointPosition = temp1.Position(6);
    T = getTransform(robot, config, 'tool0');
    positionData(i,:) = [T(1,end), T(2,end), T(3,end)];
end


ts1 = timeseries(bagselect1);
time1Data = ts1.Time-startT;


end

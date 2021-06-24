function colorOutput = cmuColor(printType)
% Returns a 3x1 vector of RBG values for the requested official CMU color 
% given in 'printType'. No argument will return 'red-web'. The currently 
% supported colors areprimary colors 'red-web', 'red-print', 'gray', and 
% 'dark-gray', and secondary colors 'gold', 'teal', 'blue', 'green', and 
% 'dark-green'. 'random' will select one of these colors at random. Nathan
% added 'carnegie-red','black','steel-gray','iron-gray','scots-rose','gold-thread','green-thread',
% 'teal-thread','blue-thread','sky-blue','hall-tan','brick-beige',
% 'hornbostel-teal','palladian-green','weaver-blue','skibo-red'

if nargin == 0
    printType = 'web';
end

colorList = {'red-web', 'red-print', 'gray', 'dark-gray', 'gold',...
    'teal', 'blue', 'green', 'dark-green','carnegie-red','black',...
    'steel-gray','iron-gray','scots-rose','gold-thread','green-thread',...
    'teal-thread','blue-thread','sky-blue','hall-tan','brick-beige',...
    'hornbostel-teal','palladian-green','weaver-blue','skibo-red'};

if strcmp(printType, 'random')
    printType = colorList{randi(length(colorList))};
end

switch printType
    case 'red-web'
        colorOutput = 1/255*[187 0 0];
    case 'red-print'
        colorOutput = 1/255*[176 28 46];
    case 'gray'
        colorOutput = 1/255*[244 244 244];
    case 'dark-gray'
        colorOutput = 1/255*[102 102 102];
    case 'gold'
%         colorOutput = 1/255*[170 102 0];
colorOutput = 1/255*[255 173 0];
    case 'teal'
        colorOutput = 1/255*[0 102 119];
    case 'blue'
        colorOutput = 1/255*[34 68 119];
    case 'green'
        colorOutput = 1/255*[18 220 0];
%         colorOutput = 1/255*[0 136 85];
    case 'dark-green'
        colorOutput = 1/255*[34 68 51];
    %%% Nathan's colors
    % Main colors
    case 'carnegie-red'
        cmyk_val = [0 100 79 20];
        colorOutput = cmyk(cmyk_val);
    case 'black'
        cmyk_val = [0 0 0 100];
        colorOutput = cmyk(cmyk_val);
    case 'steel-gray'
        cmyk_val = [0 0 0 30];
        colorOutput = cmyk(cmyk_val);
    case 'iron-gray'
        cmyk_val = [0 0 0 70];
        colorOutput = cmyk(cmyk_val);
    % Secondary colors
    case 'scots-rose'
        cmyk_val = [0 92 72 0];
        colorOutput = cmyk(cmyk_val);
    case 'gold-thread'
        cmyk_val = [0 32 100 0];
        colorOutput = cmyk(cmyk_val);
    case 'green-thread'
        cmyk_val = [92 2 100 12];
        colorOutput = cmyk(cmyk_val);
    case 'teal-thread'
        cmyk_val = [100 0 40 20];
        colorOutput = cmyk(cmyk_val);
    case 'blue-thread'
        cmyk_val = [100 80 6 32];
        colorOutput = cmyk(cmyk_val);
    case 'sky-blue'
        cmyk_val = [100 10 3 16];
        colorOutput = cmyk(cmyk_val);
    % campus palette
    case 'hall-tan'
        cmyk_val = [15 15 30 15];
        colorOutput = cmyk(cmyk_val);
    case 'brick-beige'
        cmyk_val = [10 11 23 0];
        colorOutput = cmyk(cmyk_val);
    case 'hornbostel-teal'
        cmyk_val = [85 50 58 41];
        colorOutput = cmyk(cmyk_val);
    case 'palladian-green'
        cmyk_val = [60 25 45 0];
        colorOutput = cmyk(cmyk_val);
    case 'weaver-blue'
        cmyk_val = [97 84 44 40];
        colorOutput = cmyk(cmyk_val);
    case 'skibo-red'
        cmyk_val = [15 100 87 35];
        colorOutput = cmyk(cmyk_val);
end

 
        function rgb_val = cmyk(cmyk_val)
            cc = cmyk_val(1)/100; mm = cmyk_val(2)/100; 
            yy = cmyk_val(3)/100; kk = cmyk_val(4)/100; 
            R = (1-cc)*(1-kk);
            G = (1-mm)*(1-kk);
            B = (1-yy)*(1-kk);
            rgb_val = [R,G,B];
        end
end
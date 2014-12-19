function [dataMat] = TrackPoints(filepath, filenames)
%CinePlex Data Analysis Function
%Parses text files from Cineplex to determine time and angle
%Plots each tracked point connected to the others with lines
%Tracked points should be in columns in the order that you would like
%those points connected in the figures

%Inputs: filepath (string) excluding file name
%        filenames (cell array of strings) file names separated by commas
%        numVars (int) number of points tracked

%Meghan Jimenez
%SINAPSE
%8 December 2014

%Clear anything in the figure
clf

%Global variable to track frame number through slider callback
global f

for c = 1:length(filenames)
    
%Construct file name
filename = cell2mat(strcat(filepath,'/',filenames(c),'.txt'));

%scans the data file for comma delimited strings
dataMat = load(filename);

%Calculates the number of variables (columns) and data points (rows)
[Length, numVars] = size(dataMat);

%Initialize matrix for holding theta values
thetaMat(1:Length, 1:(numVars - 2)/2 - 2) = 0;

%Set mins and maxes to the first of each x and y
xmin = dataMat(1,3);
xmax = dataMat(1,3);
ymin = dataMat(2,3);
ymax = dataMat(2,3);

%Find max and min for x and y by checking through each tracked point's 
%x and y values and comparing
for j = 3:numVars
    cmax = max(dataMat(:,j));
    cmin = min(dataMat(:,j));
    
    if mod(j,2) == 0
        if cmax > ymax
            ymax = cmax;
        end
        if cmin < ymin
            ymin = cmin;
        end
    
    else
        if cmax > ymax
            xmax = cmax;
        end
        if cmin < ymin
            xmin = cmin;
        end
    end      
        
end

%Loop over all of the variables in each frame of the video
%Skips time and frame columns, looks at variables in x,y pairs
%Tracks the previous (ppx & ppy, px & py) and current (x, y) point to
%find the angle between them and create a theta matrix
lineLen = 0;
px = 0;
py = 0;
for i = 1:Length
    for j = 3:numVars - 2
        if mod(j,2) ~= 0
            ppx = px; %2-back x,y
            ppy = py;
            
            px = dataMat(i,j); %Previous x, y
            py = dataMat(i,j + 1);

            x = dataMat(i,j + 2); %Current x, y
            y = dataMat(i,j + 3);
            
            %Store the previous line length
            pLineLen = lineLen;
            
            %Calculate the distance between points
            dx = x - px;
            dy = y - py;
            lineLen = sqrt(dx^2 + dy^2);
            
            %Calculate the angle with law of cosines
            if j > 3
                hypDx = ppx - x;
                hypDy = ppy - y;
                hyp = sqrt(hypDx^2 + hypDy^2);
                theta = acosd((pLineLen^2 + lineLen^2 - hyp^2)...
                   /(2*pLineLen*lineLen));
                thetaMat(i,(j-1)/2 - 1) = theta;
            end            
        end
    end
end

length(thetaMat)

end

%Plotting and display
%Setup color ordering
ColOrd = get(gca,'ColorOrder');
%Initialize variables for tracking frame (f) and characters pressed (CH)
f = 1;
CH = 0;
% setup GUI (slider)
hFig = figure();
%hAx = axes('Parent',hFig);
uicontrol('Parent',hFig, 'Style','slider', 'Value',1, 'Min',1,...
    'Max',Length, 'SliderStep',[1 1], ...
    'Position',[150 5 300 20], 'Callback',@slider_callback); 
%Define axis
%axis([xmin,xmax,ymin,ymax]);

%Plot lines to show angle between joints
while f <= Length
    for j = 3:numVars - 2
        if mod(j,2) ~= 0
            hold on
            
            px = dataMat(f,j); %Previous x, y
            py = dataMat(f,j + 1);

            x = dataMat(f,j + 2); %Current x, y
            y = dataMat(f,j + 3);
            
            %Plot lines
            line([px, x],[py, y], 'Color', ColOrd(j,:));            
            title(strcat('Stimulated Rat Limb Motion @ t = ', ...
                num2str( dataMat(f,2)), ', f = ', num2str(dataMat(f,1))))
            if j > 3
                text(px, py + (ymax - ymin)/25, ...
                strcat(num2str(thetaMat(f,(j-1)/2 - 1)), '°'))
            end
            
            hold off

        end
        
    end
                
    
    %Navigation through plots
    %Currently you have to click after moving the slider to update plot
    if j > 3 %This probably needs to be fixed for multi-version
        
        w = waitforbuttonpress();

        if w
            CH = get(gcf,'CurrentCharacter');
        end
        
        if CH == 98 %b navigates back
            if f ~= 1
                f = f - 1;
            end
        end
        if CH == 110 %n navigates forward
            if f < Length
                f = f + 1;
            end
        end
        if CH == 32 %space breaks
            break
        end
        if CH == 102 %f allows text input to navigate by frame number
            frame = input('Frame:');
            if frame > Length
                frame = Length;
            end
            f = frame;
        end
        
        if CH == 116 %t allows text input to navigate by time
            time = input('Time (s):');
            tFrame = find(abs((dataMat(:,2) - time)) < 0.017);
            f = tFrame(1);
        end
        
        %clears axis to draw the next figure
        cla 
        
    end
end

%clears figure
clf

%plots angle between points over time
plot(dataMat(:,2), thetaMat)
title('Rat Leg Angular Motion')
xlabel('Time (s)')
ylabel('Angle (degrees)')

end

%Callback function for slider
function sVal = slider_callback(hObj, ~ ) 
global f %global frame count must be initialized

%get the slider value
sVal = round(get(hObj,'Value'));
%If the frame number isn't the same as the slider value, change the frame
%number
if sVal ~= f
    f = sVal;
end
end
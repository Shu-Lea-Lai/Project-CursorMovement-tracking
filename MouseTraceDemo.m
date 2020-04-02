    whichScreen = Screen('OpenWindow',0,[175 175 175],[100 100 900 400]);%max(Screen('Screens'));
    [win,rect] = Screen('OpenWindow', whichScreen, [175 175 175],[100 100 900 400]);

    % Move the cursor to the center of the screen
    theX = round(rect(RectRight) / 2);
    theY = round(rect(RectBottom) / 2);
    SetMouse(theX,theY,whichScreen);
    ShowCursor;

    % Wait for a click and hide the cursor
%     Screen(win,'FillRect',0);
    Screen(win,'DrawText','Drag mouse (i.e. hold button down) to draw',50,50,255);
    Screen('Flip', win);
    while (1)
        [x,y,buttons] = GetMouse(win);
        if buttons(1) || KbCheck
          break;
        end
    end
    Screen(win,'DrawText','Drag mouse (i.e. hold button down) to draw',50,50,0);

    % Loop and track the mouse, drawing the contour
    [theX,theY] = GetMouse(win);
    thePoints = [theX theY];
    Screen(win,'DrawLine',255,theX,theY,theX,theY);
    Screen('Flip', win);
    sampleTime = 0.01;
    startTime = GetSecs;
    nextTime = startTime+sampleTime;
    while 1
        [x,y,buttons] = GetMouse(win);
        if ~buttons(1)
            break;
        end
        if (x ~= theX || y ~= theY)
            [numPoints, two]=size(thePoints);
            for i= 1:numPoints-1
                Screen(win,'DrawLine',128,thePoints(i,1),thePoints(i,2),thePoints(i+1,1),thePoints(i+1,2));
            end
            Screen('Flip', win);
            theX = x; theY = y;
        end
        if (GetSecs > nextTime)
            thePoints = [thePoints ; x y];
            nextTime = nextTime+sampleTime;
        end
    end
clear all;
debug = 2; % mode 2: developer mode
mylaptop = 0;
Screen('Preference', 'SkipSyncTests', 1);

try
    if mylaptop==1
        Screen('Preference', 'SkipSyncTests', 1);
    end
    if (debug==1|debug==2)
        [whichScreen, rect] = Screen('OpenWindow',0,[175 175 175],[100 100 900 400]);%max(Screen('Screens'));
        
    else %if i'm testing
        [win rect] = Screen('OpenWindow',0,[175 175 175]);
    end

    %% instructions
    textsize=30;
    rest_index='Rest Index Fingers on F and J';
    Screen('DrawText',whichScreen,rest_index,50,50,50); 
    Screen('Flip',whichScreen);
    
    legal=0;
    if debug ~= 0
        while legal == 0
            [keydown secs keycode]=KbCheck;
            key=KbName(keycode);
            if strcmp(key,'space')
                legal=1; 
            end
        end
    elseif debug == 1
        legal=1;
    end
    

    %% Move the cursor to the center of the screen
    [win,rect2] = Screen('OpenWindow', whichScreen, [175 175 175],[100 100 900 400]); 
    theX = round(rect(RectRight) / 2);
    theY = round(rect(RectBottom) / 2);
    SetMouse(theX,theY,whichScreen);
    ShowCursor;

    %% Wait for a click and hide the cursor
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
    
    
    
    %% clear screen 
    Screen('DrawText',win,"endnow",50,50,50);
    Screen('Flip',win); 
    WaitSecs(1)
    sca;
    
    plot(thePoints(:,1),theRect(RectBottom)-thePoints(:,2));
catch
    
%     sca;  
end  
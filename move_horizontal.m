clear all;
debug = 2;
mylaptop = 0;
Screen('Preference', 'SkipSyncTests', 1);

image_location=[pwd '/images/']; 
move_button = [image_location 'button.jpeg'];

try
    if mylaptop==1
        Screen('Preference', 'SkipSyncTests', 1);
    end
    if (debug==1 || debug==2)
        [whichScreen, rect] = Screen('OpenWindow',0,[175 175 175],[10 10 850 400]);%max(Screen('Screens'));
        
    else %if i'm testing
        [whichScreen rect] = Screen('OpenWindow',0,[175 175 175]);
    end
    
    %% get values
    [center(1), center(2)] = RectCenter(rect);
    winWidth = rect(RectRight);
    winHeight = rect(RectBottom);
    XLowLimit = ceil(winWidth*1/5);
    XUpLimit = ceil(winWidth*4/5);
    Y_pos_bar = ceil(winHeight*6/10);
    Y_pos_dot = ceil(winHeight* (6/10) * (8.5/10));
    dotcolor = 255;
    dotsize = 10;
    %% instructions
%     textsize=30;
%     rest_index='Rest Index Fingers on F and J';
%     Screen('DrawText',whichScreen,rest_index,50,50,50); 
%     Screen('Flip',whichScreen);
%     
%     legal=0;
%     if debug ~= 0
%         while legal == 0
%             [keydown secs keycode]=KbCheck;
%             key=KbName(keycode);
%             if strcmp(key,'space')
%                 legal=1; 
%             end
%         end
%     elseif debug == 1
%         legal=1;
%     end
    

    %% Move the cursor to the center of the screen
    if (debug==1 || debug==2)
        [win,rect2] = Screen('OpenWindow', whichScreen, [175 175 175],[10 10 850 400]);
        
    else %if i'm testing
        [win,rect2] = Screen('OpenWindow',whichScreen,[175 175 175]);
    end
    
     
    theX = round(winWidth/ 2);
    theY = round(winHeight / 2);
    SetMouse(XLowLimit,Y_pos_bar,whichScreen);
    ShowCursor;

    %% Wait for a click and hide the cursor
    Screen('DrawText',win, 'Drag mouse (i.e. hold button down) to draw',round(winWidth*1/6),round(winWidth*1/3), 255);
    picture=Screen('MakeTexture',win,imread(move_button));
    Screen('DrawTexture',win,picture,[],[XLowLimit-15,Y_pos_bar-15 ,XLowLimit+15 ,Y_pos_bar+15 ]);
    Screen('DrawDots', win, [XLowLimit, Y_pos_dot], dotsize, dotcolor);
    Screen('Flip', win);
    Screen('Close',picture);
    
    while (1)
        [x,y,buttons] = GetMouse(win);
        if buttons(1) || KbCheck
          break;
        end
    end
    Screen(win,'DrawText','Start now',round(winWidth*1/6),round(winWidth*1/3), 255);
    Screen('Flip',win);
%     HideCursor(); 
        
    %% Move dot and Track Cursor
    
    
    [theX,theY] = GetMouse(win);
    theCursors = [theX theY]; % org start
%     theDots = [XLowLimit y_pos_dot];
    
    sampleTime = 0.01;
    startTime = GetSecs;
    nextTime = startTime+sampleTime;
    
    loop=0;
    
    while 1
        
        loop;
        [x,y,buttons] = GetMouse(win);
        vdot = sin(GetSecs);
        xdot_i = XLowLimit + vdot*GetSecs;
        buttons(1)
        if ~buttons(1) %if 
            break;
        end
        if (x ~= theX || y ~= theY) %if moves
            [numPoints, two]=size(theCursors);
            currentXpos = theCursors(numPoints,1);

            picture=Screen('MakeTexture',win,imread(move_button));
            if (currentXpos > XUpLimit)
                
                Screen('DrawDots', win, [XUpLimit, Y_pos_dot], dotsize, dotcolor);
                Screen('DrawTexture',win,picture,[],[XUpLimit-15,Y_pos_bar-15, XUpLimit+15 ,Y_pos_bar+15 ]);
            elseif (currentXpos < XLowLimit)
                
                Screen('DrawDots', win, [XLowLimit, Y_pos_dot], dotsize, dotcolor);
                Screen('DrawTexture',win,picture,[],[XLowLimit-15,Y_pos_bar-15 ,XLowLimit+15 ,Y_pos_bar+15 ]);
            else
                
                Screen('DrawDots', win, [xdot_i, Y_pos_dot], dotsize, dotcolor);%, center, 1)
                Screen('DrawTexture',win,picture,[],[currentXpos-15,Y_pos_bar-15 ,currentXpos+15 ,Y_pos_bar+15 ]);
            end
            
            Screen('Close',picture);
            Screen('Flip', win);
            
            theX = x; theY = y;
        end
        if (GetSecs > nextTime)
            theCursors = [theCursors ; x y];
%             theDots = [Xdot ; xdot_i y_pos_dot]
            nextTime = nextTime+sampleTime;
        end
        loop = loop + 1
    end
    
    
    
    %% clear screen 
    Screen('DrawText',win,"windnow",round(winWidth*1/6),round(winWidth*1/3), 50);
    Screen('Flip',win); 
    WaitSecs(1)
    sca;
    
    plot(theCursors(:,1),theRect(RectBottom)-theCursors(:,2));
catch
    
%     sca;  
end  
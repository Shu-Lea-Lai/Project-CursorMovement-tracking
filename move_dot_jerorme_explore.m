
%% initial values
debug =1;
mylaptop = 0;
% Screen('Preference', 'SkipSyncTests', 1);

image_location=[pwd '/images/']; 
data_location=[pwd '/data/'];
subid=input(' subject # ');
move_button = [image_location 'arrow.png'];
steady_bar = [image_location 'bar.png'];

if debug==0
    filename=[data_location 'simple' num2str(subid) '.txt'];
    matdatas = [data_location 'simple' num2str(subid)];
else
    filename=[data_location 'debug_simple' num2str(subid) '.txt'];
    matdatas =  [data_location 'debug_simple' num2str(subid)];
end

%% create window
if mylaptop==1
    Screen('Preference', 'SkipSyncTests', 1);
end
if (debug==1 || debug==2)
    [whichScreen, rect] = Screen('OpenWindow',0,[175 175 175],[10 10 850 400]);%max(Screen('Screens'));
    [win,rect2] = Screen('OpenWindow', whichScreen, [175 175 175],[10 10 850 400]);
else %if i'm testing
    [whichScreen, rect] = Screen('OpenWindow',0,[175 175 175]);
    [win,rect2] = Screen('OpenWindow',whichScreen,[175 175 175]);
end

%% Set various parameters of experiments

% values of system

[center(1), center(2)] = RectCenter(rect);
winWidth = rect(RectRight);
winHeight = rect(RectBottom);
winMiddel = round(winWidth/2);
XleftLimit = ceil(winWidth*1/20);
XrightLimit = ceil(winWidth*19/20);
Xmiddel = ceil((XleftLimit + XrightLimit)/2);
DotMoveRange = round(winWidth*7/10);
Y_pos_bar = ceil(winHeight*6/10);
Y_pos_dot = ceil(winHeight* (6/10) * (8.5/10));
dotcolor = 255;
dotsize = 10;

% values of trial
numBlock = 1;
numTrialpB = 1;
timelastpT = 5; % time (seconds) that a trial last

%% instructions

WaitSpace(win, 'Tracking object. Press space to start', debug);
WaitSecs(0.1);


%% Exp start
store_line = [];
fid=fopen(filename,'wt');
for iblock = 1:numBlock
    
    store_line_trial = [];
    
    %% parameters of the movmenet function
    
    b = 0;
    a = 0;
    c = randi([1,10],1);
    
    while a >= b || b/a>2 %we don't want when the ratio is more than 2, or when a is bigger than a
        
        b = randi([1,10],1);
        a = randi([1,10],1);
    end
    
    for itrial = 1:numTrialpB
        
        %% set a random starting position (time)  
        start_t = randi([0 40],1);
        input_t = start_t;
        tspan = [0 100]; %time from 0 to 100, garentees to have periods included 
        theta_sol = ode45(@(t,y) (b - a*sin(y)), tspan, start_t);
        Xinitial = GenXdot(theta_sol,input_t,DotMoveRange,c,winMiddel);

        
        
        %% Move mouse to the start position
        SetMouse(Xinitial,Y_pos_bar,whichScreen);
        ShowCursor;
        
        %% Wait for a click and hide the cursor
        Screen('DrawText',win, 'Click on mouse to start',round(winWidth*1/6),round(winWidth*1/3), 255);
        
        Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
        [img_bar, ~, alpha_bar] = imread(steady_bar);
        [img_button, ~, alpha_button] = imread(move_button);
        img_bar(:, :, 4) = alpha_bar;
        img_button(:,:,4) = alpha_button;
        pic_bar = Screen('MakeTexture', win, img_bar);
        pic_button = Screen('MakeTexture',win,img_button);
        
        textureflip = @(loc, pic) Screen('DrawTexture',win,pic,[],loc);
        
        textureflip([XleftLimit,Y_pos_bar-10 ,XrightLimit ,Y_pos_bar+10 ], pic_bar); %initial position
        textureflip([Xinitial-15,Y_pos_bar-15 ,Xinitial+15 ,Y_pos_bar+15 ], pic_button);
        Screen('DrawDots', win, [Xinitial, Y_pos_dot], dotsize, dotcolor);
        Screen('Flip', win);
        Screen('Close',pic_button);
        Screen('Close',pic_bar);

        while 1
            [MouseX,MouseY,buttons] = GetMouse(win);
            if buttons(1) || KbCheck
              break;
            end
        end
        Screen(win,'DrawText','Start now',round(winWidth*1/6),round(winWidth*1/3), 255);
        Screen('Flip',win);
        HideCursor ; 

        %% Move dot and Track Cursor


        [theX,theY] = GetMouse(win);
        theCursors = [theX theY]; % org start
        theDots = [Xinitial Y_pos_dot];
        theFlips = 0;
%         accTime = [];
        accTime = 0;

        flipTime = 0.001;
        startTime = GetSecs;
        nextFlip = startTime+flipTime;
        loop=1;
        innerloop = 1;  
        time0 = GetSecs;
        xdot_i = theDots(1); 
             
        pic_button = Screen('MakeTexture',win,img_button);
        pic_bar = Screen('MakeTexture', win, img_bar);
        while 1
            
            loopstart = GetSecs;
%             x = checktime(time0,1);
            x = GetSecs;
            loop;
            [MouseX, MouseY, ~] = GetMouse(win);
            
            
            if ((GetSecs - time0) >= timelastpT)
                break;
            end
            
            x = checktime(x,2);
            if (GetSecs > nextFlip)
                
%                 x = checktime(x,3);
                %% movment of x
                innerloop;
                CurrentBetweenTime = GetSecs - time0;
                input_t = CurrentBetweenTime + start_t;
                MovePerUnit = 0.005*winWidth; %constant to move
                
                xdot_i = GenXdot(theta_sol,input_t,DotMoveRange,c,winMiddel);
                
                [numPoints, two] = size(theCursors);
                currentXpos = theCursors(numPoints,1);
                 
                %% Draw slide arrow

                
                x = checktime(x,4);

                if (currentXpos > XrightLimit)
                    
                    textureflip([XleftLimit,Y_pos_bar-10 ,XrightLimit ,Y_pos_bar+10 ], pic_bar);
                    textureflip([XrightLimit-15,Y_pos_bar-15, XrightLimit+15 ,Y_pos_bar+15 ], pic_button);
                elseif (currentXpos < XleftLimit)
                    
                    textureflip([XleftLimit,Y_pos_bar-10 ,XrightLimit ,Y_pos_bar+10 ], pic_bar);
                    textureflip([XleftLimit-15,Y_pos_bar-15 ,XleftLimit+15 ,Y_pos_bar+15 ], pic_button);
                else
                    
                    textureflip([XleftLimit,Y_pos_bar-10 ,XrightLimit ,Y_pos_bar+10 ], pic_bar);
                    textureflip([currentXpos-15,Y_pos_bar-15 ,currentXpos+15 ,Y_pos_bar+15 ], pic_button);
                end
                x = checktime(x,5);
                
                Screen('DrawDots', win, [xdot_i, Y_pos_dot], dotsize, dotcolor);%, center, 1)
%                 Screen('Close',pic_button);
%                 Screen('Close',pic_bar);
                Screen('Flip', win);
                x = checktime(x,6);

                theCursors = [theCursors ; MouseX MouseY]; %#ok<AGROW>
                theDots = [theDots ; xdot_i Y_pos_dot]; %#ok<AGROW>
                theFlips = [theFlips; round(GetSecs - loopstart,4)];%#ok<AGROW>
                accTime = [accTime; round(GetSecs - startTime,4)];%#ok<AGROW>
                nextFlip = nextFlip + flipTime;
                innerloop = innerloop+1;
                
                x = checktime(x,7);
            end
            x = checktime(x,8);
            loop = loop + 1;
            [~] = checktime(loopstart,"end");
        end
        
        %% store data into store_line
        len_pos = max(size(theCursors,1), size(theDots,1));
        if size(theCursors,1) ~= size(theDots,1) 
            fprintf("%s","unmatched length of pos");
            system('kill');
        end
        
        curtime = split(string(datetime('now'))," ");
        repq = @(x) repmat(x, len_pos, 1);
        currentline = [repq(curtime(1)) repq(curtime(2)) repq(iblock), repq(itrial), theCursors, theDots];
        if iblock == 1 && itrial == 1
            store_line_trial = ["Day" "Time" "Block" "Trial" "CursorX" "CursorY" "DotX" "DotY"];
        end
        store_line_trial = [store_line_trial; currentline]; %#ok<AGROW>
        
        %% print data into .txt
        for isubline = 1:len_pos 
            fprintf(fid,'%12s', curtime(1), curtime(2));
            fprintf(fid,'%5d',iblock, itrial, a, b, c);
            fprintf(fid,'%.4f', theFlips(isubline));
            fprintf(fid,'%2s'," ");
            fprintf(fid,'%.4f', accTime(isubline));
            fprintf(fid,'%10d', theCursors(isubline, 1), theCursors(isubline, 2),...
                theDots(isubline, 1), theDots(isubline, 2)); 
            fprintf(fid,'\n');
        end
    end %endtrial
    
    %% print inter-Block message
    if debug ~= 1 && iblock ~= numBlock 
        WaitSpace(win,'Block is ended, hiting space to continue', debug);
        fprintf('something went wrong\n')
    elseif debug ~= 1  && iblock == numBlock
        WaitSpace(win, 'You reach the end of the experiment, informing the experimenter before leaving', debug);
    end
    WaitSecs(0.2);
    store_line = [store_line; store_line_trial]; %#ok<AGROW>
    
    
end %endblock

fclose(fid);
save(matdatas);
Screen('CloseAll');

%% Functions
function WaitSpace(win, text, debug)

%     textsize=30;
    Screen('DrawText',win,text,50,50,50); 
    Screen('Flip',win);

    legal=0;
    if debug ~= 1 % hit space to continue 
        while legal == 0
            [~, ~, keycode]=KbCheck;
            key = KbName(keycode);
            if strcmp(key,'space')
                legal=1; 
            end
        end
    end
end

function [xdot_i] = GenXdot(theta_sol,input_t,DotMoveRange,c,winMiddel)

    width_factor = round(DotMoveRange/2);
    theta_t = deval(theta_sol, input_t);
    r_t = width_factor * cos(c*input_t);
    x_t = r_t* cos(theta_t); %range:[-r_t,r_t]
    xdot_i = winMiddel + x_t;   
    xdot_i = round(xdot_i);
  
end

function x = checktime(time0,point)

    between = GetSecs - time0;
    x = GetSecs;
    if between*1000>10
        [point,between*1000];
    end
end

    
    
    
    
    
    
    
    
    
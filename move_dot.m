
clear all;
%% initial values
debug = 0;
mylaptop = 0;
Screen('Preference', 'SkipSyncTests', 1);

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

%% get values

% values of system

[center(1), center(2)] = RectCenter(rect);
winWidth = rect(RectRight);
winHeight = rect(RectBottom);
XLowLimit = ceil(winWidth*1/5);
XUpLimit = ceil(winWidth*4/5);
Y_pos_bar = ceil(winHeight*6/10);
Y_pos_dot = ceil(winHeight* (6/10) * (8.5/10));
dotcolor = 255;
dotsize = 10;

% values of trial
numBlock = 1;
numTrialpB = 1;
timelastpT = 20;

%% instructions

WaitSpace(win, 'Tracking object. Press space to start', debug);
WaitSecs(0.1);


%% Exp start
store_line = [];
fid=fopen(filename,'wt');
for iblock = 1:numBlock
    
    store_line_trial = [];
    
    for itrial = 1:numTrialpB

        %% Move mouse to the start position
        SetMouse(XLowLimit,Y_pos_bar,whichScreen);
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
        
        textureflip([XLowLimit,Y_pos_bar-10 ,XUpLimit+15 ,Y_pos_bar+10 ], pic_bar);
        textureflip([XLowLimit-15,Y_pos_bar-15 ,XLowLimit+15 ,Y_pos_bar+15 ], pic_button);
        Screen('DrawDots', win, [XLowLimit, Y_pos_dot], dotsize, dotcolor);
        Screen('Flip', win);
        Screen('Close',pic_button);
        Screen('Close',pic_bar);

        while 1
            [x,y,buttons] = GetMouse(win);
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
        theDots = [XLowLimit Y_pos_dot];

        sampleTime = 0.01;
        flipTime = 0.01;
        startTime = GetSecs;
        nextTime = startTime+sampleTime;
        nextFlip = startTime+flipTime;
        loop=1;
        time0 = GetSecs*1000;
        xdot_i = theDots(1);

        toRight = 1; 
        while 1

            loop;
            [x,y,buttons] = GetMouse(win);

            if ((GetSecs - time0/1000) >= timelastpT)
                break;
            end

            if (GetSecs > nextFlip)

                timei = GetSecs*1000;
                timeBetween = timei - time0;
                
                
                
                vdot = ceil((sin(timeBetween/100 + (rand-0.5)*4 )+0.5) * 20); %x0

%                 make sure vdot doesn't pass boundary if vdot is negative
                if ((xdot_i <= XLowLimit && vdot < 0) || (xdot_i >= XUpLimit && vdot < 0))
                    while (vdot<0)
                        vdot = ceil((sin(timeBetween/100 + (rand-0.5)*4 )+0.5) * 10);
                    end
                end

%                 Bouncing back if reach boundary
                if (xdot_i <= XLowLimit )
                    toRight = 1;                
                elseif (xdot_i >= XUpLimit)
                    toRight = 0;                 
                end

                if toRight ==1 
                    xdot_i = xdot_i + vdot;
                else
                    xdot_i = xdot_i - vdot;
                end

                [numPoints, two]=size(theCursors);
                currentXpos = theCursors(numPoints,1);
                 
                % Draw slide arrow
                pic_button = Screen('MakeTexture',win,img_button);
                pic_bar = Screen('MakeTexture', win, img_bar);
                
                

                if (currentXpos > XUpLimit)
                    textureflip([XLowLimit,Y_pos_bar-10 ,XUpLimit+15 ,Y_pos_bar+10 ], pic_bar);
                    textureflip([XUpLimit-15,Y_pos_bar-15, XUpLimit+15 ,Y_pos_bar+15 ], pic_button);
                elseif (currentXpos < XLowLimit)
                    textureflip([XLowLimit,Y_pos_bar-10 ,XUpLimit+15 ,Y_pos_bar+10 ], pic_bar);
                    textureflip([XLowLimit-15,Y_pos_bar-15 ,XLowLimit+15 ,Y_pos_bar+15 ], pic_button);
                else
                    textureflip([XLowLimit,Y_pos_bar-10 ,XUpLimit+15 ,Y_pos_bar+10 ], pic_bar);
                    textureflip([currentXpos-15,Y_pos_bar-15 ,currentXpos+15 ,Y_pos_bar+15 ], pic_button);
                end

                Screen('DrawDots', win, [xdot_i, Y_pos_dot], dotsize, dotcolor);%, center, 1)
                Screen('Close',pic_button);
                Screen('Close',pic_bar);
                Screen('Flip', win);

                theX = x; theY = y;
                theCursors = [theCursors ; x y]; %#ok<AGROW>
                theDots = [theDots ; xdot_i Y_pos_dot]; %#ok<AGROW>
                nextFlip = nextFlip + flipTime;
            end
            loop = loop + 1;
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
            fprintf(fid,'%5d',iblock, itrial);
            fprintf(fid,'%10d', theCursors(isubline, 1), theCursors(isubline, 2), theDots(isubline, 1), theDots(isubline, 2));            
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

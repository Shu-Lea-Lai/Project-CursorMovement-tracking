% PsychDefaultSetup(2); 
Screen('Preference', 'SkipSyncTests', 1);
% Screen('Preference', 'VisualDebugLevel', 3);
[win, rect] = Screen('OpenWindow', 0, [100 100 100], [0 0 400 400]);

Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[img_bar, ~, alpha] = imread('bar.png');
img_bar(:, :, 4) = alpha;
pic_bar = Screen('MakeTexture', win, img_bar);
Screen('DrawTexture', win, pic_bar, [], [0 0 250 250]);
Screen('Flip', win);
WaitSecs(1);
sca;
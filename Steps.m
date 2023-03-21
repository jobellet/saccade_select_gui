%% 
%%1) change this according to your needs
file_name = 'subject1';
file_x = ['C:\Users\joach\Dropbox\saccade_select_GUI/',file_name,'.mat']; % Path of the file containing the horizontal eye position

file_y = file_x; % Path of the file containing the vertical eye position

horizontal_field = 'eye.x';%'myfieldname';% name of the variable in the file_x.mat file that contains the horizontal eye position

vertical_field = 'eye.y';%'myotherfieldname'; % name of the variable in the file_.mat file that contains the horizontal eye position

save_path = ['C:\Users\joach\Dropbox\saccade_select_GUI/',file_name,'_labels.mat']; % Path of the file that will contain the labels of saccades


%% 2) Launch the GUI
saccade_select_gui('horizontal_file_path',file_x,...
    'horizontal_field',horizontal_field,...
'vertical_file_path',file_y,...
'vertical_field',vertical_field,... 
'save_path',save_path,...
'ylim',[-1 1])


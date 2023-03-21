function varargout = saccade_select_gui(varargin)
% List of name value argument
%     'horizontal_file_path'  <--- Path of the file containing the horizontal eye position
%     'vertical_file_path'  <--- Path of the file containing the vertical eye position
%     'label_file_path' <--- Path of the file containing saccade labels
%
%     'horizontal_field'
%     'vertical_field'
%     'label_field'
%     'Pred_prob_file_path'
%     'Pred_prob_field'
%     'save_path'
%     'ylim'


% SACCADE_SELECT_GUI MATLAB code for saccade_select_gui.fig
%      SACCADE_SELECT_GUI, by itself, creates a new SACCADE_SELECT_GUI or raises the existing
%      singleton*.
%
%      H = SACCADE_SELECT_GUI returns the handle to a new SACCADE_SELECT_GUI or the handle to
%      the existing singleton*.
%
%      SACCADE_SELECT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SACCADE_SELECT_GUI.M with the given input arguments.
%
%      SACCADE_SELECT_GUI('Property','Value',...) creates a new SACCADE_SELECT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before saccade_select_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to saccade_select_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help saccade_select_gui

% Last Modified by GUIDE v2.5 08-Aug-2020 15:52:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @saccade_select_gui_OpeningFcn, ...
    'gui_OutputFcn',  @saccade_select_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before saccade_select_gui is made visible.
function saccade_select_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to saccade_select_gui (see VARARGIN)

% Choose default command line output for saccade_select_gui
handles.output = hObject;
handles.zoom=0;
%%%%%%%%%%%%

handles.colors=[0,0,0;
    60, 180, 75;
    255, 225, 25;
    0, 130, 200;
    245, 130, 48;
    145, 30, 180;]./255;
handles.global_view=1;
handles.samples_num=500;
handles.samples_num_precision=50;
handles.min_event_size=1;
handles.min_gap_size=1;
handles.label_num=1;
handles.threshold=2;
handles.savepath=[];
handles.current_labels=ones(1,6);
handles.row=1;
handles.interv=1:handles.samples_num;
handles.interv_duration=handles.interv(end)-handles.interv(1)+1;
handles.ylim = [-1 1];
handles.undo = [];
% go through gui inputs
label_found=0;
feature_found=0;
if strcmp(varargin{1}(end-2:end),'edf') || strcmp(varargin{1}(end-2:end),'EDF') %option 1: load directly an EDF file
    handles.savepath=[varargin{1}(1:end-4),'_labels']; % <------ modify this if you want to change the name of the file saved
    handles.trial_type='cell';
    Trials= edfImport(varargin{1}(1:end-4), [1 1 1], ''); % this is the function from Alexander Pastukhov you need to have Eyelink API installed to make it work
    
    for i=1:length(varargin) % loop through function inputs to extract crucial informations
        if strcmp(varargin{i},'screen_resolution')
            handles.screen_resolution=varargin{i+1};
        end
        if strcmp(varargin{i},'ylim')
            handles.ylim=varargin{i+1};
        end
    end
    
    
    handles.X=cell(length(Trials),1);
    handles.Y=cell(length(Trials),1);
    handles.explored_interv=cell(length(Trials),1);
    handles.velocity=cell(length(Trials),1);
    handles.labels_manual=cell(length(Trials),1);
    handles.saccades_auto=cell(length(Trials),1);
    for i = 1 : length(Trials)
        handles.X{i}=Trials(i).Samples.gx(2,:);
        handles.Y{i}=Trials(i).Samples.gy(2,:);
        handles.X{i}=handles.X{i}-handles.screen_resolution(1)/2;
        handles.Y{i}=handles.Y{i}-handles.screen_resolution(2)/2;
        handles.X{i}=handles.X{i}./Trials(i).Samples.rx;
        handles.Y{i}=handles.Y{i}./Trials(i).Samples.ry;
        handles.explored_interv{i}=false(size(handles.Y{i}));
        if  label_found==0
            handles.labels_manual{i}=zeros(size(handles.Y{i}));
        else
            handles.labels_manual{i}=round(handles.labels_manual{i});
        end
        if feature_found==0
            handles.velocity{i}=[0,sqrt(diff(handles.X{i}).^2+diff(handles.Y{i}).^2)];
        end
        handles.saccades_auto{i}=false(size(handles.Y{i}));
    end
    
    
else
    for i=1:length(varargin) % loop through function inputs to extract
        if strcmp(varargin{i},'horizontal_file_path')
            if any(cellfun(@(x) ~isempty(strfind(x,'horizontal_field')),varargin))
                handles.X=load_x_y_or_label(varargin{i+1},varargin{find(cellfun(@(x) ~isempty(strfind(x,'horizontal_field')),varargin))+1});
            else
                handles.X=load_x_y_or_label(varargin{i+1},[]);
            end
            
        end
        if strcmp(varargin{i},'vertical_file_path')
            if any(cellfun(@(x) ~isempty(strfind(x,'vertical_field')),varargin))
                handles.Y=load_x_y_or_label(varargin{i+1},varargin{find(cellfun(@(x) ~isempty(strfind(x,'vertical_field')),varargin))+1});
            else
                handles.Y=load_x_y_or_label(varargin{i+1},[]);
            end
        end
    end
    
    if iscell(handles.X)
        handles.trial_type='cell';
        
        % rotate cell array if not in the adequate dimension
        if size(handles.Y,1)<size(handles.Y,2)
            handles.Y=handles.Y';
        end
        if size(handles.X,1)<size(handles.X,2)
            handles.X=handles.X';
        end
    else
        handles.trial_type='matrix';
    end
    for i=1:length(varargin) % loop through function inputs
        
        if strcmp(varargin{i},'label_file_path')
            label_found=1;
            if any(cellfun(@(x) ~isempty(strfind(x,'label_field')),varargin))
                handles.labels_manual=load_x_y_or_label(varargin{i+1},varargin{find(cellfun(@(x) ~isempty(strfind(x,'label_field')),varargin))+1});
            else
                handles.labels_manual=load_x_y_or_label(varargin{i+1},[]);
            end
            for row = 1:size(handles.Y,1)
                handles=cluster_of_labels(handles,row);
            end
            
        end
        if strcmp(varargin{i},'Pred_prob_file_path')
            feature_found=1;
            if any(cellfun(@(x) ~isempty(strfind(x,'Pred_prob_field')),varargin))
                handles.velocity=load_x_y_or_label(varargin{i+1},varargin{find(cellfun(@(x) ~isempty(strfind(x,'Pred_prob_field')),varargin))+1});
            else
                handles.velocity=load_x_y_or_label(varargin{i+1},[]);
                
            end
        end
        if strcmp(varargin{i},'save_path')
            handles.savepath=varargin{i+1};
        end
        if strcmp(varargin{i},'ylim')
            handles.ylim=varargin{i+1};
        end
    end
    if iscell(handles.X)
        handles.explored_interv=cell(size(handles.X,1),1);
        for i =1 :size(handles.X,1)
            
            if  label_found==0
                handles.labels_manual{i}=zeros(size(handles.Y{i}));
            else
                handles.labels_manual{i}=round(handles.labels_manual{i});
            end
            if feature_found==0
                handles.velocity{i}=[0,sqrt(diff(handles.X{i}).^2+diff(handles.Y{i}).^2)];
            end
        end
        
    else
        
        handles.X=double(handles.X);
        handles.Y=double(handles.Y);
        handles.explored_interv=false(size(handles.Y));
        if  label_found==0
            handles.labels_manual=zeros(size(handles.Y));
        else
            handles.labels_manual=round(handles.labels_manual);
        end
        if feature_found==0;
            handles.velocity=[zeros(size(handles.X,1),1),sqrt(diff(handles.X,1,2).^2+diff(handles.Y,1,2).^2)];
        end
        handles.saccades_auto=false(size(handles.Y));
    end
    
end



handles=cluster_of_labels(handles,handles.row);
% This sets up the initial plot

plot_everything(handles)

% set samples string
guidata(hObject, handles);


% UIWAIT makes saccade_select_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = saccade_select_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Next.
function Next_Callback(hObject, eventdata, handles)
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles; handles = release_undo_memory(handles);
if handles.global_view
    handles.interv=handles.interv+round((handles.interv(end)-handles.interv(1))/2);
    if strcmp(handles.trial_type,'cell')
        trial_length=length( handles.X{handles.row});
    else
        trial_length=size(handles.X,2);
    end
    if handles.interv(end)>trial_length
        if  handles.interv(end)==trial_length+round((handles.interv(end)-handles.interv(1))/2)
            
            if handles.row<size(handles.X,1)
                handles.row=handles.row+1;
                handles=cluster_of_labels(handles,handles.row);
                handles.interv=1:length(handles.interv);
            else
                handles.interv=(trial_length-length(handles.interv)+1):trial_length;
                disp('YOU HAVE REACHED THE END OF THE FILE')
            end
        else
            handles.interv=(trial_length-length(handles.interv)+1):trial_length;
        end
    end
else % if we are in precision mode, go to the next event
    
    handles.current_labels(handles.label_num)=handles.current_labels(handles.label_num)+1;
    if sum(cellfun(@(x) length(x),handles.clusters{handles.label_num})) < handles.current_labels(handles.label_num)
        handles.current_labels(handles.label_num)=sum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
        disp('YOU HAVE REACHED THE LAST EVENT')
    end
end

if strcmp(handles.trial_type,'cell')
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    handles.explored_interv(handles.row,handles.interv)=true;
end

plot_everything(handles)

guidata(hObject, handles);

% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
if handles.global_view
    handles.interv=handles.interv-round((handles.interv(end)-handles.interv(1))/2);
    if handles.interv(1)<1
        if  handles.interv(1)==1-round((handles.interv(end)-handles.interv(1))/2)
            if handles.row>1
                handles.row=handles.row-1;
                handles=cluster_of_labels(handles,handles.row);
                if strcmp(handles.trial_type,'cell')
                    trial_length=length( handles.X{handles.row});
                else
                    trial_length=size(handles.X,2);
                end
                handles.interv=(trial_length-length(handles.interv)+1):trial_length;
            end
        else
            handles.interv=1:length(handles.interv);
        end
    end
else
    handles.current_labels(handles.label_num)=handles.current_labels(handles.label_num)-1;
    if handles.current_labels(handles.label_num)<1
        handles.current_labels(handles.label_num)=1;
    end
end
if strcmp(handles.trial_type,'cell')
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    handles.explored_interv(handles.row,handles.interv)=true;
end
plot_everything(handles)
guidata(hObject, handles);


% --- Executes on button press in SaveLabels.
function SaveLabels_Callback(hObject, eventdata, handles)
% hObject    handle to SaveLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Wait don''t close the GUI!')
if ~isempty(handles.savepath)
    label=handles.labels_manual;
    try % this repace the variable if it already exists in the .mat file
        save(handles.savepath,'label','-append'); %'YOUR_VARIABLE_NAME'
    catch % if the variable doesn't exist yet it is going to be created
        save(handles.savepath,'label'); %'YOUR_VARIABLE_NAME'
    end
else
    labels=handles.labels_manual;
    save('sacade_select_gui_output','labels')
end
disp('SAVED')
disp('You can close the GUI now')
% --- Executes on button press in addlabel.
function addlabel_Callback(hObject, eventdata, handles)
% hObject    handle to addlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
axes(handles.axes1);
[x,~]=ginput(1);
x=round(x);
if handles.global_view==1 % here the button is used to add an event
    row=handles.row;
    interv=handles.interv;
    if x<max(interv)
        if strcmp(handles.trial_type,'cell')
            trial_length=length( handles.X{handles.row});
            handles.labels_manual{row}(max(x-10,1):min(x+10,trial_length))=handles.label_num;
            handles=cluster_of_labels(handles,row);
            handles.explored_interv{handles.row}(handles.interv)=true;
        else
            trial_length=size(handles.X,2);
            handles.labels_manual(row,max(x-10,1):min(x+10,trial_length))=handles.label_num;
            handles=cluster_of_labels(handles,row);
            handles.explored_interv(handles.row,handles.interv)=true;
        end
        
        plot_everything(handles)
        guidata(hObject, handles);
    end
else  % here the button is used to fine-tune saccade start
    all_events=cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
    row = find(all_events>=handles.current_labels(handles.label_num),1,'first');
    if row>1
        clust_end=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)}(end);
    else
        clust_end= handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)}(end);
    end
    if x<clust_end
        %delete previous label
        if row>1
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)};
        else
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)};
        end
        if strcmp(handles.trial_type,'cell')
            handles.labels_manual{row}(interv)=0;
        else
            handles.labels_manual(row,interv)=0;
        end
        % replace with modified onset
        if row>1
            handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)}=x:clust_end;
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)};
        else
            handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)}=x:clust_end;
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)};
        end
        if strcmp(handles.trial_type,'cell')
            handles.labels_manual{row}(interv)=handles.label_num;
        else
            handles.labels_manual(row,interv)=handles.label_num;
        end
        handles=cluster_of_labels(handles,row);
        if strcmp(handles.trial_type,'cell')
            handles.explored_interv{handles.row}(handles.interv)=true;
        else
            handles.explored_interv(handles.row,handles.interv)=true;
        end
        plot_everything(handles)
        guidata(hObject, handles);
    end
end



% --- Executes on button press in removelabel.
function removelabel_Callback(hObject, eventdata, handles)
% hObject    handle to removelabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
axes(handles.axes1);
[x,~]=ginput(1);
x=round(x);
if handles.global_view==1 % here the button is used to remove an event
    handles=remove_labels(hObject,handles,x);
else  % here the button is used to fine-tune event offset
    all_events=cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
    row = find(all_events>=handles.current_labels(handles.label_num),1,'first');
    if row>1
        clust_start=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)}(1);
    else
        clust_start= handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)}(1);
    end
    if x>clust_start
        %delete previous label
        if row>1
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)};
        else
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)};
        end
        if strcmp(handles.trial_type,'cell')
            handles.labels_manual{row}(interv)=0;
        else
            handles.labels_manual(row,interv)=0;
        end
        % replace with modified onset
        if row>1
            handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)}=clust_start:x;
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)};
        else
            handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)}=clust_start:x;
            interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)};
        end
        if strcmp(handles.trial_type,'cell')
            handles.labels_manual{row}(interv)=handles.label_num;
        else
            handles.labels_manual(row,interv)=handles.label_num;
        end
        handles=cluster_of_labels(handles,row);
        if strcmp(handles.trial_type,'cell')
            handles.explored_interv{handles.row}(handles.interv)=true;
        else
            handles.explored_interv(handles.row,handles.interv)=true;
        end
        plot_everything(handles)
        guidata(hObject, handles);
    end
end

% --- Executes on button press in DetectSacAuto.
function DetectSacAuto_Callback(hObject, eventdata, handles)
% hObject    handle to DetectSacAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
if strcmp(handles.trial_type,'cell')
    
    all_velocities=cat(2,handles.velocity{:});
    for i=1: length(handles.X)
        handles.saccades_auto{i}=handles.velocity{i}>handles.threshold*median(all_velocities(all_velocities~=0))/0.6745;
        handles.labels_manual{i}(find(handles.saccades_auto{i}))=1;
    end
    
    temp=handles.label_num;
    handles.label_num=1;
    for row = 1:length(handles.Y)
        handles = merge_labels(handles,row);
        handles = cluster_of_labels(handles,row);
    end
    handles.label_num=temp;
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    handles.labels_manual(:,:)=0;
    all_velocities=handles.velocity(:);
    handles.saccades_auto=handles.velocity>handles.threshold*median(all_velocities(all_velocities~=0))/0.6745;
    handles.labels_manual(find(handles.saccades_auto))=1;
    
    temp=handles.label_num;
    handles.label_num=1;
    for row = 1:size(handles.Y,1)
        handles=merge_labels(handles,row);
        handles=cluster_of_labels(handles,row);
    end
    handles.label_num=temp;
    handles.explored_interv(handles.row,handles.interv)=true;
end
plot_everything(handles)
guidata(hObject, handles)

function labelnumber_Callback(hObject, eventdata, handles)
% hObject    handle to labelnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labelnumber as text
%        str2double(get(hObject,'String')) returns contents of labelnumber as a double
handles.undo = handles;
handles.label_num=str2double(get(hObject,'String'));
if handles.label_num<1
    handles.label_num=1;
end
plot_everything(handles)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function labelnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_everything(handles)

if handles.global_view==0
    all_events=cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
    row = find(all_events>=handles.current_labels(handles.label_num),1,'first'); % the trial is changed if no event is found in the current trial
    if row>1
        interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)};
    else
        interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)};
    end
    if strcmp(handles.trial_type,'cell')
        interv=max(interv(1)-handles.samples_num_precision,1):min(interv(end)+handles.samples_num_precision,length(handles.Y{row}));
    else
        interv=max(interv(1)-handles.samples_num_precision,1):min(interv(end)+handles.samples_num_precision,size(handles.Y,2));
    end
else
    row=handles.row;
    interv=handles.interv;
end
if strcmp(handles.trial_type,'cell')
    temp_labels=unique(handles.labels_manual{row});
else
    temp_labels=unique(handles.labels_manual(row,:));
end

%%%%%% plot horizontal eye position
axes(handles.axes1);
if strcmp(handles.trial_type,'cell')
    plot(interv,handles.X{row}(interv))
else
    plot(interv,handles.X(row,interv))
end
hold on
for i=1:length(temp_labels)
    if temp_labels(i)~=0 %%% draw the dots on event
        if strcmp(handles.trial_type,'cell')
            label_index=find(handles.labels_manual{row}==temp_labels(i));
            label_index=label_index(ismember(label_index,interv));
            if handles.label_num==temp_labels(i)
                plot(label_index,handles.X{row}(label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','MarkerSize',10,'LineStyle','none')
            else
                plot(label_index,handles.X{row}(label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','LineStyle','none')
            end
        else
            label_index=find(handles.labels_manual(row,:)==temp_labels(i));
            label_index=label_index(ismember(label_index,interv));
            if handles.label_num==temp_labels(i)
                plot(label_index,handles.X(row,label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','MarkerSize',10,'LineStyle','none')
            else
                plot(label_index,handles.X(row,label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','LineStyle','none')
            end
        end
    end
end
hold off
set(gca,'xticklabel',{[]})
if handles.zoom==1
    ylim('auto');
else
    if ~any(isnan(handles.ylim)) %control zoom of Y axis
        ylim(handles.ylim)
    end
end

xlim([min(interv),max(interv)])


%%%%%% plot vertical eye position
axes(handles.axes2);
if strcmp(handles.trial_type,'cell')
    plot(interv,handles.Y{row}(interv))
else
    plot(interv,handles.Y(row,interv))
end
hold on
for i=1:length(temp_labels)
    if temp_labels(i)~=0 %%% draw the dots on event
        if strcmp(handles.trial_type,'cell')
            label_index=find(handles.labels_manual{row}==temp_labels(i));
            label_index=label_index(ismember(label_index,interv));
            if handles.label_num==temp_labels(i)
                plot(label_index,handles.Y{row}(label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','MarkerSize',10,'LineStyle','none')
            else
                plot(label_index,handles.Y{row}(label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','LineStyle','none')
            end
        else
            label_index=find(handles.labels_manual(row,:)==temp_labels(i));
            label_index=label_index(ismember(label_index,interv));
            if handles.label_num==temp_labels(i)
                plot(label_index,handles.Y(row,label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','MarkerSize',10,'LineStyle','none')
            else
                plot(label_index,handles.Y(row,label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','LineStyle','none')
            end
        end
    end
end
hold off
set(gca,'xticklabel',{[]})
if handles.zoom==1
    ylim('auto');
else
    if ~any(isnan(handles.ylim)) %control zoom of Y axis
        ylim(handles.ylim)
    end
end


xlim([min(interv),max(interv)])


%%%%%%% Plot velocity
axes(handles.axes3);

if strcmp(handles.trial_type,'cell')
    plot(interv,handles.velocity{row}(interv))
else
    plot(interv,handles.velocity(row,interv))
end

hold on
for i=1:length(temp_labels)
    if temp_labels(i)~=0
        if strcmp(handles.trial_type,'cell')
            label_index=find(handles.labels_manual{row}==temp_labels(i));
            label_index=label_index(ismember(label_index,interv));
            if handles.label_num==temp_labels(i)
                plot(label_index,handles.velocity{row}(label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','MarkerSize',10,'LineStyle','none')
            else
                plot(label_index,handles.velocity{row}(label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','LineStyle','none')
            end
        else
            label_index=find(handles.labels_manual(row,:)==temp_labels(i));
            label_index=label_index(ismember(label_index,interv));
            
            if handles.label_num==temp_labels(i)
                plot(label_index,handles.velocity(row,label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','MarkerSize',10,'LineStyle','none')
            else
                plot(label_index,handles.velocity(row,label_index),'color',handles.colors(temp_labels(i),:),'Marker','.','LineStyle','none')
            end
        end
    end
end
hold off
if handles.zoom==1
    ylim('auto');
else
end
xlim([min(interv),max(interv)])

%%%%% Plot explored trials
axes(handles.axes12);
if strcmp(handles.trial_type,'cell')
    temp=zeros(length(handles.explored_interv),100);
    for i=1:length(handles.explored_interv)
        temp_len=length(handles.explored_interv{i});
        temp_len=round(find(handles.explored_interv{i})*100/temp_len);
        temp_len=unique(temp_len(temp_len>0 & temp_len<101));
        temp(i,temp_len)=1;
    end
    imagesc(temp);
    colormap bone
    hold on
    plot(round(handles.interv(1)*100/length(handles.explored_interv{i})),handles.row,'b.')
else
    imagesc(handles.explored_interv);
    colormap bone
    hold on
    plot(handles.interv(1),handles.row,'b.')
end
hold off
ax = gca;
% ax.Visible = 'off';
ylabel('Trial number')
xlabel('Time sample')
% update window range
set(handles.trial_nb,'String',num2str(handles.row))
set(handles.sample_pos,'String',[num2str(handles.interv(1)),'-',num2str(handles.interv(end))])

function handles=merge_labels(handles,row)
if strcmp(handles.trial_type,'cell')
    clusters=bwconncomp(handles.labels_manual{row}~=handles.label_num);
else
    clusters=bwconncomp(handles.labels_manual(row,:)~=handles.label_num);
end
clusters=clusters.PixelIdxList;
clusters=clusters(cellfun(@(x) length(x)<handles.min_gap_size, clusters));
if strcmp(handles.trial_type,'cell')
    for i=1:length(clusters)
        handles.labels_manual{row}(clusters{i})=handles.label_num;
    end
else
    for i=1:length(clusters)
        handles.labels_manual(row,clusters{i})=handles.label_num;
    end
end
function handles=cluster_of_labels(handles,row)
% group continuous samples of a same label into a cluster
if strcmp(handles.trial_type,'cell')
    temp=unique(handles.labels_manual{row});
else
    temp=unique(handles.labels_manual(row,:));
end
temp=temp(temp~=0);
for j=1:length(temp)
    if strcmp(handles.trial_type,'cell')
        clusters=bwconncomp(handles.labels_manual{row}==temp(j));
    else
        clusters=bwconncomp(handles.labels_manual(row,:)==temp(j));
    end
    clusters=clusters.PixelIdxList;
    handles.clusters{temp(j)}{row}=clusters(cellfun(@(x) length(x)>handles.min_event_size, clusters));
    if strcmp(handles.trial_type,'cell')
        handles.labels_manual{row}(handles.labels_manual{row}==temp(j))=0;
    else
        handles.labels_manual(row,handles.labels_manual(row,:)==temp(j))=0;
    end
    for i=1:length(handles.clusters{temp(j)}{row})
        if strcmp(handles.trial_type,'cell')
            handles.labels_manual{row}(handles.clusters{temp(j)}{row}{i})=temp(j);
        else
            handles.labels_manual(row,handles.clusters{temp(j)}{row}{i})=temp(j);
        end
    end
end


function min_interv_Callback(hObject, eventdata, handles)
% hObject    handle to min_sac_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_sac_length as text
%        str2double(get(hObject,'String')) returns contents of min_sac_length as a double
handles.undo = handles;
handles.min_gap_size=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function min_interv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_interv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_sac_length_Callback(hObject, eventdata, handles)
% hObject    handle to min_sac_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_sac_length as text
%        str2double(get(hObject,'String')) returns contents of min_sac_length as a double
handles.undo = handles;
handles.min_event_size=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function min_sac_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_sac_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text4.
function text4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if handles.global_view==1
    switch eventdata.Key
        case 'u'
            undo_Callback(hObject, eventdata, handles)
        case 'a'
            addlabel_Callback(hObject, eventdata, handles)
        case 'r'
            removelabel_Callback(hObject, eventdata, handles)
        case 'space'
            precisionmode_Callback(hObject, eventdata, handles)
        case 'x'
            previous_Callback(hObject, eventdata, handles)
        case 'v'
            Next_Callback(hObject, eventdata, handles)
        case 'rightarrow'
            Next_Callback(hObject, eventdata, handles)
        case 'leftarrow'
            previous_Callback(hObject, eventdata, handles)
        case 'i'
            CheckUnexplored_Callback(hObject, eventdata, handles)
        case 'z'
            zoombutton_Callback(hObject, eventdata, handles)
        case 'g'
            startend_Callback(hObject, eventdata, handles)
    end
else % if we are in precision mode
    switch eventdata.Key
        case 'u'
            undo_Callback(hObject, eventdata, handles)
        case 's'
            addlabel_Callback(hObject, eventdata, handles) %this is to change the onset of the event
        case 'f'
            removelabel_Callback(hObject, eventdata, handles) %this is to change the offset of the event
        case 'r'
            CheckUnexplored_Callback(hObject, eventdata, handles) %this is to delete event
        case 'space'
            precisionmode_Callback(hObject, eventdata, handles)
        case 'x'
            previous_Callback(hObject, eventdata, handles)
        case 'v'
            Next_Callback(hObject, eventdata, handles)
        case 'rightarrow'
            Next_Callback(hObject, eventdata, handles)
        case 'leftarrow'
            previous_Callback(hObject, eventdata, handles)
        case 'z'
            zoombutton_Callback(hObject, eventdata, handles)
        case 'g'
            startend_Callback(hObject, eventdata, handles)
    end
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
handles.undo = handles;
handles.threshold=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in precisionmode.
function precisionmode_Callback(hObject, eventdata, handles)
% hObject    handle to precisionmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
if handles.global_view
    if isfield(handles,'clusters')
        if ~isempty(sum(cellfun(@(x) length(x),handles.clusters{handles.label_num})))
            all_events=cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
            row = handles.row; % the event is in the same trial than the interval
            if all_events(row)>0
                handles.global_view=0;
                
                % set the event viewed as the event occuring the closest to the center of the current
                % interval
                
                distance_to_center=cellfun(@(x) min(abs(x-handles.interv(round(length(handles.interv)/2)))),handles.clusters{handles.label_num}{row});
                clust_nb=find(distance_to_center==min(distance_to_center),1,'first');
                if row>1 % in case it is not the first trial
                    handles.current_labels(handles.label_num)=clust_nb+all_events(row-1); % add the number of event in previous trials
                else
                    handles.current_labels(handles.label_num)=clust_nb;
                end
                
                
                % changes the text on buttons
                set(handles.precisionmode, 'String', 'Switch to global view');
                set(handles.Next, 'String', 'Next event>>>');
                set(handles.previous, 'String', '<<<Previous event');
                set(handles.removelabel,'String','Select offset')
                set(handles.addlabel,'String','Select onset')
                set(handles.CheckUnexplored, 'String', 'Remove');
                set(handles.text12,'String','Shortcut: s')
                set(handles.text13,'String','Shortcut: f')
                set(handles.text24,'String','Shortcut: r')
            else
                disp('NO EVENT SELECTED IN THIS TRIAL')
            end
        end
    end
else
    handles.global_view=1;
    
    % set interval such that it is centered on the last event observed
    all_events=cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
    row = find(all_events>=handles.current_labels(handles.label_num),1,'first');
    handles.row=row;
    handles=cluster_of_labels(handles,handles.row);
    if row>1
        clust_start = handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)}(1);
        clust_end = handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)}(end);
    else
        clust_start = handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)}(1);
        clust_end = handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)}(end);
    end
    
    
    
    handles.interv=round(clust_start+(clust_end-clust_start)/2)+(1:handles.interv_duration)-round(handles.interv_duration/2);
    if handles.interv(1)<1
        handles.interv= handles.interv-handles.interv(1)+1;
    end
    if strcmp(handles.trial_type,'cell')
        if handles.interv(end)>length(handles.Y{row})
            handles.interv = length(handles.Y{row})-handles.interv_duration+1:length(handles.Y{row});
        end
        handles.interv=max(1,handles.interv(1)):min(length(handles.Y{row}),handles.interv(end));
    else
        if handles.interv(end)>size(handles.Y,2)
            handles.interv = size(handles.Y,2)-handles.interv_duration+1:size(handles.Y,2);
        end
        handles.interv=max(1,handles.interv(1)):min(size(handles.Y,2),handles.interv(end));
    end
    % changes the text on buttons
    set(handles.precisionmode, 'String', 'Switch to event view');
    set(handles.Next, 'String', 'Next interval>>>');
    set(handles.previous, 'String', '<<<Previous interval');
    set(handles.removelabel,'String','Remove lablel')
    set(handles.CheckUnexplored, 'String', 'Sample of interest');
    set(handles.addlabel,'String','Add label')
    set(handles.text12,'String','Shortcut: a')
    set(handles.text13,'String','Shortcut: r')
    set(handles.text24,'String','Shortcut: i')
end
if strcmp(handles.trial_type,'cell')
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    handles.explored_interv(handles.row,handles.interv)=true;
end
plot_everything(handles)
guidata(hObject, handles);


function M=load_x_y_or_label(filename,field)
% - write the full path of filename which corresponds to the
% horizontal eye position, the vertical eye position or the label of the eye
% movement. It can be a .csv or a .mat file.

% -varargin is the field for accessing the variable in .mat file
% file. For instance: 'data.eye_position.horizontal'
% Don't need to specify it if there is only one one matrice in the .mat file

if strcmp(filename(end-2:end),'csv')
    M =csvread(filename);
elseif strcmp(filename(end-2:end),'mat')
    M =load(filename);
    variables=fields(M);
    if ~isempty(field)
        fields_names=field;
        dot_pos=strfind(fields_names,'.');
        for i=1:length(dot_pos)
            M=M.(fields_names(1:dot_pos(1)-1));
            fields_names=fields_names(dot_pos(1)+1:end);
            dot_pos=strfind(fields_names,'.');
        end
        M=M.(fields_names);
    else
        M=M.(variables{1});
    end
end

if size(M,2)==1
    M=M';
end

function handles=remove_labels(hObject,handles,x)
handles.undo = handles;
if ~isempty(x) % if in global view
    row=handles.row;
    nearest_clust=NaN(1,length(handles.clusters));
    distance=NaN(1,length(handles.clusters));
    for i=1:length(handles.clusters)
        if ~isempty(handles.clusters{i})
            if  length(handles.clusters{i})>=row
                if ~isempty(handles.clusters{i}{row})
                    min_dist=cellfun(@(n) min(abs(n-x)),handles.clusters{i}{row});
                    nearest_clust(i)=find(min_dist==min(min_dist));
                    distance(i)=min(min_dist);
                end
            end
        end
        
    end
    if exist('min_dist','var')
        
        if strcmp(handles.trial_type,'cell')
            handles.labels_manual{row}(handles.clusters{find(distance==min(distance),1,'first')}{row}{nearest_clust(find(distance==min(distance),1,'first'))})=0;
            handles=cluster_of_labels(handles,row);
            handles.explored_interv{handles.row}(handles.interv)=true;
        else
            handles.labels_manual(row,handles.clusters{find(distance==min(distance),1,'first')}{row}{nearest_clust(find(distance==min(distance),1,'first'))})=0;
            handles=cluster_of_labels(handles,row);
            handles.explored_interv(handles.row,handles.interv)=true;
        end
        plot_everything(handles)
    end
    
else % if in precision mode
    all_events=cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
    row = find(all_events>=handles.current_labels(handles.label_num),1,'first');
    if row>1
        interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)-all_events(row-1)};
    else
        interv=handles.clusters{handles.label_num}{row}{handles.current_labels(handles.label_num)};
    end
    if strcmp(handles.trial_type,'cell')
        handles.labels_manual{row}(interv)=0;
        handles.explored_interv{row}(interv)=true;
    else
        handles.labels_manual(row,interv)=0;
        handles.explored_interv(row,interv)=true;
    end
    plot_everything(handles)
    pause(.5)
    handles=cluster_of_labels(handles,row);
    plot_everything(handles)
end
if handles.current_labels(handles.label_num)>cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}))
    handles.current_labels(handles.label_num)=cumsum(cellfun(@(x) length(x),handles.clusters{handles.label_num}));
end
guidata(hObject, handles);



function smoothingedit_Callback(hObject, eventdata, handles)
% hObject    handle to smoothingedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothingedit as text
%        str2double(get(hObject,'String')) returns contents of smoothingedit as a double
handles.undo = handles;
new_val=get(hObject,'String');
if ~strcmp(new_val,'NaN')
    new_val = round(str2double(new_val));
    if mod(new_val,2)==0
        new_val=new_val+1;
    end
    if new_val<2
        new_val=3;
    end
    set(handles.smoothingedit,'String',num2str(new_val));
    set(hObject,'Value',new_val)
    disp(new_val)
    if strcmp(handles.trial_type,'cell')
        for i = 1 : length(handles.X)
            handles.velocity{i}=[0,sqrt(diff(sgolayfilt(handles.X{i}',2,new_val)',1,2).^2+diff(sgolayfilt(handles.Y{i}',2,new_val)',1,2).^2)];
        end
    else
        handles.velocity=[zeros(size(handles.X,1),1),sqrt(diff(sgolayfilt(handles.X',2,new_val)',1,2).^2+diff(sgolayfilt(handles.Y',2,new_val)',1,2).^2)];
    end
else
    handles.velocity=[zeros(size(handles.X,1),1),sqrt(diff(handles.X,1,2).^2+diff(handles.Y,1,2).^2)];
end
if strcmp(handles.trial_type,'cell')
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    handles.explored_interv(handles.row,handles.interv)=true;
end
plot_everything(handles)
guidata(hObject, handles);
set(hObject, 'Enable', 'off');
drawnow;
set(hObject, 'Enable', 'on');

% --- Executes during object creation, after setting all properties.
function smoothingedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothingedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CheckUnexplored.
function CheckUnexplored_Callback(hObject, eventdata, handles)
% hObject    handle to CheckUnexplored (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.global_view==1 % in global view this button is used for selecting an interesting interval
    labeled=handles.labels_manual>0;
    if any(~handles.explored_interv(:)&~labeled(:))
        [I,J] = ind2sub(size(handles.X),find(handles.velocity==max(handles.velocity(~handles.explored_interv&~labeled))&~handles.explored_interv&~labeled,1,'first'));
        handles.row=I;
        disp(handles.row)
        handles=cluster_of_labels(handles,handles.row);
        handles.interv=J+(-300:300);
        handles.interv=max(handles.interv(1),1):min(handles.interv(end),trial_length);
        
        handles.explored_interv(I,handles.interv)=true;
        plot_everything(handles)
        guidata(hObject, handles);
    else
        disp('Every samples have been analysed already')
    end
else % in precision mode this button is used to remove the label
    handles=remove_labels(hObject,handles,[]);
end

function trial_nb_Callback(hObject, eventdata, handles)
% hObject    handle to trial_nb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trial_nb as text
%        str2double(get(hObject,'String')) returns contents of trial_nb as a double
handles.undo = handles;
value=str2double(get(hObject,'String'));
if value<=size(handles.Y,1)
    handles.row=value;
end
set(handles.trial_nb,'String',num2str(handles.row))
if strcmp(handles.trial_type,'cell')
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    handles.explored_interv(handles.row,handles.interv)=true;
end
handles=cluster_of_labels(handles,handles.row);
plot_everything(handles)
guidata(hObject, handles);
set(hObject, 'Enable', 'off');
drawnow;
set(hObject, 'Enable', 'on');

% --- Executes during object creation, after setting all properties.
function trial_nb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trial_nb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sample_pos_Callback(hObject, eventdata, handles)
% hObject    handle to sample_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_pos as text
%        str2double(get(hObject,'String')) returns contents of sample_pos as a double
handles.undo = handles;
val=get(hObject,'String');
first_sample=str2double(val(1:strfind(val,'-')-1));
last_sample=str2double(val(strfind(val,'-')+1:end));

first_sample=max(first_sample,1);

if strcmp(handles.trial_type,'cell')
    mxlength=length(handles.Y{handles.row});
    if first_sample>length(handles.Y{handles.row})
        first_sample = max(mxlength-300,1);
    end
else
    mxlength=size(handles.Y,2);
    if first_sample>size(handles.Y,2)
        first_sample=max(mxlength-300,1);
    end
end

handles.interv_duration=min(last_sample,mxlength)-first_sample+1;
handles.interv=first_sample:min(last_sample,mxlength);
if strcmp(handles.trial_type,'cell')
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    handles.explored_interv(handles.row,handles.interv)=true;
end
plot_everything(handles)
guidata(hObject, handles);
set(hObject, 'Enable', 'off');
drawnow;
set(hObject, 'Enable', 'on');

% --- Executes during object creation, after setting all properties.
function sample_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ID1.
function ID1_Callback(hObject, eventdata, handles)
% hObject    handle to ID1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
handles.label_num=1;
plot_everything(handles)
guidata(hObject, handles);

% --- Executes on button press in ID2.
function ID2_Callback(hObject, eventdata, handles)
% hObject    handle to ID2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
handles.label_num=2;
plot_everything(handles)
guidata(hObject, handles);


% --- Executes on button press in ID3.
function ID3_Callback(hObject, eventdata, handles)
% hObject    handle to ID3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
handles.label_num=3;
plot_everything(handles)
guidata(hObject, handles);


% --- Executes on button press in ID4.
function ID4_Callback(hObject, eventdata, handles)
% hObject    handle to ID4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
handles.label_num=4;
plot_everything(handles)
guidata(hObject, handles);


% --- Executes on button press in ID5.
function ID5_Callback(hObject, eventdata, handles)
% hObject    handle to ID5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
handles.label_num=5;
plot_everything(handles)
guidata(hObject, handles);


% --- Executes on button press in ID6.
function ID6_Callback(hObject, eventdata, handles)
% hObject    handle to ID6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
handles.label_num=6;
plot_everything(handles)
guidata(hObject, handles);


% --- Executes on button press in zoombutton.
function zoombutton_Callback(hObject, eventdata, handles)
% hObject    handle to zoombutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
handles.zoom = mod(handles.zoom+1,2);
guidata(hObject, handles);
plot_everything(handles)


% --- Executes on button press in startend.
function startend_Callback(hObject, eventdata, handles)
% hObject    handle to startend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.undo = handles;
[x,~]=ginput(2);
x=sort(x);
x=round(x);
row= handles.row;
if strcmp(handles.trial_type,'cell')
    trial_length=length( handles.X{handles.row});
    handles.labels_manual{row}(max(x(1),1):min(x(2),trial_length))=handles.label_num;
    handles=cluster_of_labels(handles,row);
    handles.explored_interv{handles.row}(handles.interv)=true;
else
    trial_length=size(handles.X,2);
    handles.labels_manual(row,max(x(1),1):min(x(2),trial_length))=handles.label_num;
    handles=cluster_of_labels(handles,row);
    handles.explored_interv(handles.row,handles.interv)=true;
end

plot_everything(handles)
guidata(hObject, handles);


% --- Executes on button press in undo.
function undo_Callback(hObject, eventdata, handles)
% hObject    handle to undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.undo)
    handles = handles.undo;
    plot_everything(handles)
end
guidata(hObject, handles);

function handles = release_undo_memory(handles)
if isfield(handles,'undo')
    if isfield(handles.undo,'undo')
        if isfield(handles.undo.undo,'undo')
            handles.undo.undo.undo = [];
        end
    end
end



% --- Executes on button press in Ready4UnEye.
function Ready4UnEye_Callback(hObject, eventdata, handles)
% hObject    handle to Ready4UnEye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.trial_type,'cell')
    try
        X2save = cell2mat(handles.X);
    catch
        X2save = cell2mat(handles.X');
    end
    
    try
        Y2save = cell2mat(handles.Y);
    catch
        Y2save = cell2mat(handles.Y');
    end
    
    try
        L2save = cell2mat(handles.labels_manual);
    catch
        L2save = cell2mat(handles.labels_manual');
    end
else
    X2save = handles.X;
    Y2save = handles.Y;
    L2save = handles.labels_manual;
    
    
end
X2save = transpose_when_needed(X2save);
Y2save = transpose_when_needed(Y2save);
L2save = transpose_when_needed(L2save);


selpath = uigetdir;
if isfile([selpath,'\X.csv'])
    answer = questdlg('What do you want to do?', ...
        'Training files already exist', ...
        'Combine them','Replace them','Cancel','Cancel');
    % Handle response
    switch answer
        case 'Combine them'
            tempX = csvread([selpath,'\X.csv']);
            tempY = csvread([selpath,'\Y.csv']);
            tempL = csvread([selpath,'\L.csv']);
            tempX = transpose_when_needed(tempX);
            tempY = transpose_when_needed(tempY);
            tempL = transpose_when_needed(tempL);
            
            answer2 = 'Combine them anyway';
            if ~isempty(strfind(tempX(1,:),round(X2save(1,:),5,'significant')))
                answer2 = questdlg('What do you want to do?', ...
                    'Same data is already present in the original file', ...
                    'Combine them anyway','Cancel','Cancel');
            end
            if strcmp( answer2,'Combine them anyway')
                if size(tempX,1) == 1
                    
                    X2save = [tempX,X2save];
                else
                    X2save = [tempX;X2save];
                end
                
                if size(tempY,1) == 1
                    Y2save = [tempY,Y2save];
                else
                    Y2save = [tempY;Y2save];
                end
                
                if size(tempL,1) == 1
                    L2save = [tempL,L2save];
                else
                    L2save = [tempL;L2save];
                end
                
                Ready4UnEye(selpath,X2save,Y2save,L2save)
            end
        case 'Replace them'
            Ready4UnEye(selpath,X2save,Y2save,L2save)
        case 'Cancel'
            disp('Training file was not created')
    end
else
    Ready4UnEye(selpath,X2save,Y2save,L2save)
end

function Ready4UnEye(selpath,X2save,Y2save,L2save)
    csvwrite([selpath,'\X.csv'],X2save)
    csvwrite([selpath,'\Y.csv'],Y2save)
    csvwrite([selpath,'\L.csv'],L2save)
    
function array = transpose_when_needed(array)
    if size(array,2) == 1
        array = array';
    end
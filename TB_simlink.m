
function varargout = responsespectrum(varargin)
% RESPONSESPECTRUM MATLAB code for responsespectrum.fig
%      RESPONSESPECTRUM, by itself, creates a new RESPONSESPECTRUM or raises the existing
%      singleton*.
%
%      H = RESPONSESPECTRUM returns the handle to a new RESPONSESPECTRUM or the handle to
%      the existing singleton*.
%
%      RESPONSESPECTRUM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESPONSESPECTRUM.M with the given input arguments.
%
%      RESPONSESPECTRUM('Property','Value',...) creates a new RESPONSESPECTRUM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before responsespectrum_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to responsespectrum_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help responsespectrum
% Last Modified by GUIDE v2.5 16-Feb-2016 12:50:25
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @responsespectrum_OpeningFcn, ...
                   'gui_OutputFcn',  @responsespectrum_OutputFcn, ...
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
% --- Executes just before responsespectrum is made visible.
function responsespectrum_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to responsespectrum (see VARARGIN)
% Choose default command line output for responsespectrum
handles.output = hObject;
clc
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes responsespectrum wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = responsespectrum_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
% --- Executes during object creation, after setting all properties.
function unitms2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unitms2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in openfile.
function openfile_Callback(hObject, eventdata, handles)
filechoice=get(handles.openfile,'Value');
    switch filechoice
        case 2 
            [FileName,Path]=uigetfile({'*.*'},'Open File');
            if isequal (FileName, 0)
                errordlg('Invalid file type','Error message');
                return
            else
                FileName=fullfile(Path,FileName);
                fid=fopen(FileName,'r');
                THM=fscanf(fid,'%g %g', [2 inf]);
                fclose(fid);
                THM=THM';
            end    
        case 3
            THM=inputdlg('Enter the matrix name','File Preloaded into Matlab');
    end
 setappdata(handles.openfile,'THM',THM);
 t=double(THM(:,1));
 y=double(THM(:,2));
    
    tmx=max(t);
    tmi=min(t);
    n=length(y);
    dt=(tmx-tmi)/(n-1);
    sr=1/dt;
setappdata(handles.openfile,'sr',sr);
    
    %nsamples
    N=strcat(num2str(n),' samples');
    set(handles.nsamples,'String',N);
    guidata(hObject, handles);
    
    % SR
    A=strcat('SR=',num2str(sr),' samples/sec');
    set(handles.SR,'String',A);
    guidata(hObject, handles);
    
    % dt
    B=strcat('dt=',num2str(dt),' sec');
    set(handles.dt,'String',B);
    guidata(hObject, handles);
    guidata(hObject,handles)
        
% hObject    handle to openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns openfile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from openfile
% --- Executes during object creation, after setting all properties.
function openfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in algsel.
function algsel_Callback(hObject, eventdata, handles)
algorithm=get(handles.algsel,'Value');
switch algorithm
    case 2  %Kelly Richman
        ialgorithm=1; 
    case 3  %Smallwood
        ialgorithm=2; 
    case 4  %Newmark
        ialgorithm=3;
end
%load variables (data)
THM = getappdata(handles.openfile , 'THM');
fn = getappdata(handles.freq0, 'fn');
Q = getappdata(handles.Qsol , 'Q');
damp= getappdata(handles.dampsol, 'damp');
iunit = getappdata(handles.unitms2 , 'iunit');
horlabelt = getappdata(handles.HlabelT , 'horlabelt');
 t=double(THM(:,1));
 y=double(THM(:,2));
    
    tmx=max(t);
    tmi=min(t);
    n=length(y);
    dt=(tmx-tmi)/(n-1);
    sr=1/dt;
    u0=0; % data for newmark
    udot0=0;
  	gam = 1/2; 
    beta = 1/4;
    
% (copied from SRS.m openfile)
tmax=(tmx-tmi) + 1/fn(1);
limit = round( tmax/dt);
n=limit;
yy=y;  % define the variable, no values in it yet
%  SRS & (ADDITIONAL) Newmark engine
for j=1:1000
    omega=2.*pi*fn(j);
    omegad=omega*sqrt(1.-(damp^2));
    cosd=cos(omegad*dt);
    sind=sin(omegad*dt);
    domegadt=damp*omega*dt;
    if ialgorithm==3 % Newmark
            omegad = 2*damp*omega;
            %c = 2*damp*omega*m;
            c = 2*damp*omega; % c divided by m
            n=length(y);
            wn=2*pi*fn;
            %p=-m*y; %instead of using p we will use y, the acceleration
            %given in the openfile we load
            % kgor = k + gam/(beta*dt)*c + m/(beta*dt^2); 
            kgor = omega^2 + gam/(beta*dt)*c + 1/(beta*dt^2); %divided by m
            % a = m/(beta*dt) + gam*c/beta; 
            a = 1/(beta*dt) + gam*c/beta;  %divided by m
            % b = 0.5*m/beta + dt*(0.5*gam/beta - 1)*c; 
            b = 0.5*1/beta + dt*(0.5*gam/beta - 1)*c; %divided by m
            dy = diff(y); % we use y (the acceleration)
            u = zeros(n,1);
            udot = u;
            u2dot = u;
            u(1) = u0;
            udot(1) = udot0;
            % u2dot(1) = 1/m*(p(1)-k*u0-c*udot0); 
            u2dot(1) = (y(1)-omega*u0-c*udot0); % since we previously divided by m
            for i = 1:(n-1)  % need to review this algorithm
                deltaA = dy(i) + a*udot(i) + b*u2dot(i);
                du_i = deltaA/kgor;
                dudot_i = gam/(beta*dt)*du_i - gam/beta*udot(i) + dt*(1-0.5*gam/beta)*u2dot(i);
                du2dot_i = 1/(beta*dt^2)*du_i - 1/(beta*dt)*udot(i) - 0.5/beta*u2dot(i);
                u(i+1) = du_i + u(i);
                udot(i+1) = dudot_i + udot(i);
                u2dot(i+1) = du2dot_i + u2dot(i);
            end
            resp=u2dot;
    else
        if ialgorithm==1 % Kelly Richman
            a1(j)=2.*exp(-domegadt)*cosd;
            a2(j)=-exp(-2.*domegadt);
            b1(j)=2.*domegadt;
            b2(j)=omega*dt*exp(-domegadt);
            b2(j)=b2(j)*( (omega/omegad)*(1.-2.*(damp^2))*sind -2.*damp*cosd );
            b3(j)=0;
        else   %Smallwood
            E=exp(-damp*omega*dt);
            K=omegad*dt;
            C=E*cos(K);
            S=E*sin(K);
            Sp=S/K;
    %
            a1(j)=2*C;
            a2(j)=-E^2;
            b1(j)=1.-Sp;
            b2(j)=2.*(Sp-C);
            b3(j)=E^2-Sp;
        end
        
        forward=[ b1(j),  b2(j),  b3(j) ];    
        back   =[     1, -a1(j), -a2(j) ];    
    %    
        resp=filter(forward,back,yy);
    end
    % we take the max and min response value for that frequency fn(j)
        x_pos(j)= max(resp);
        x_neg(j)= min(resp);
    %   
        jnum=j; 
    if  fn(j) > sr/8.
        break
    end
    fn(j+1)=fn(1)*(2. ^ (j*(1./12.)));    
end
srs_max = max(x_pos);
if max( abs(x_neg) ) > srs_max
    srs_max = max( abs(x_neg ));
end
srs_min = min(x_pos);
if min( abs(x_neg) ) < srs_min
    srs_min = min( abs(x_neg ));
end  
absx_neg=abs(x_neg);
% store solution
setappdata(handles.algsel, 'fn', fn);
setappdata(handles.algsel, 'x_pos', x_pos);
setappdata(handles.algsel, 'x_neg', x_neg);
% represent
if horlabelt==1
    if iunit==0
        handles.ejex=1./fn;
        x_neg=x_neg./9.81; 
        x_pos=x_pos./9.81; 
        plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration (G)');
    else
        handles.ejex=1./fn;
        plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration a_g(m/s^2)');
    end
else
    if iunit==0
      handles.ejex=fn;
      x_neg=x_neg./9.81; 
      x_pos=x_pos./9.81; 
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel ('Peak Acceleration (G)');
    else
      handles.ejex=fn;
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel ('Peak Acceleration a_g(m/s^2)');
    end      
end
title(sprintf('Acceleration Shock Response Spectrum Q=%g',Q));
grid on;
%set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','log');
legend ('positive','negative',2,'Location','best');
guidata(hObject,handles);
% hObject    handle to algsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns algsel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algsel
% --- Executes during object creation, after setting all properties.
function algsel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in damping.
function damping_Callback(hObject, eventdata, handles)
idamp= get(handles.damping,'Value');
switch idamp
    case 2
        damp=inputdlg('Enter damping ratio (typically 0.05)','Damping Ratio');
        damp=str2double(damp);
        Q=1/(2*damp);
    case 3
        Q=inputdlg('Enter the amplification factor (typically Q=10)','Amplification Factor Q');
        Q=str2double(Q);
        damp=1/(2*Q);
end
% stores
setappdata(handles.dampsol, 'damp', damp);
setappdata(handles.Qsol, 'Q', Q);
  % Damp
    C=strcat('damp= ',num2str(damp));
    set(handles.dampsol,'String',C);
  % Amplification Factor
    D=strcat('Q= ',num2str(Q)); 
    set(handles.Qsol,'String',D);
guidata(hObject,handles);
% hObject    handle to damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns damping contents as cell array
%        contents{get(hObject,'Value')} returns selected item from damping
% --- Executes during object creation, after setting all properties.
function damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function freq0_Callback(hObject, eventdata, handles) % starting frequency
fn(1)=str2double(get(hObject,'String'));
sr = getappdata(handles.openfile , 'sr');
if fn(1)>sr/30;
    fn(1)=sr/30; 
end
handles.fn=fn(1);
setappdata(handles.freq0, 'fn', fn)
guidata(hObject,handles);
% hObject    handle to freq0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of freq0 as text
%        str2double(get(hObject,'String')) returns contents of freq0 as a double
% --- Executes during object creation, after setting all properties.
function freq0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in psvelocity.
function psvelocity_Callback(hObject, eventdata, handles) 
fn= getappdata(handles.algsel , 'fn');
x_pos= getappdata(handles.algsel , 'x_pos');
x_neg= getappdata(handles.algsel , 'x_neg');
iunit= getappdata(handles.unitms2 , 'iunit');
horlabelt = getappdata(handles.HlabelT , 'horlabelt');
damp = getappdata(handles.dampsol , 'damp');
jnum=length(fn);
for j=1:jnum  
    if iunit==0   
       x_pos(j)=386.*x_pos(j)/(2.*pi*fn(j));
       x_neg(j)=386.*x_neg(j)/(2.*pi*fn(j));   
    else
       x_pos(j)=x_pos(j)/(2.*pi*fn(j));
       x_neg(j)=x_neg(j)/(2.*pi*fn(j));   
    end
end 
srs_max = max(x_pos);
if max( abs(x_neg) ) > srs_max
    srs_max = max( abs(x_neg ));
end
srs_min = min(x_pos);
if min( abs(x_neg) ) < srs_min
    srs_min = min( abs(x_neg ));
end  
% represent:
if horlabelt==1
    if iunit==0
        handles.ejex=1./fn;
        plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
         ylabel('Velocity (in/sec)');
    else
        handles.ejex=1./fn;
        plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel('Velocity (m/sec)'); 
    end
else
    if iunit==0
      handles.ejex=fn;
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
       ylabel('Velocity (in/sec)');
    else
      handles.ejex=fn;
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel('Velocity (m/sec)'); 
    end      
end
 
Q=1./(2.*damp);
out5 = sprintf(' Pseudo Velocity Shock Response Spectrum Q=%g ',Q);
title(out5);
grid;
%set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','log');
legend ('positive','negative',2, 'Location', 'best');
guidata(hObject,handles)
% hObject    handle to psvelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in outputfile.
function outputfile_Callback(hObject, eventdata, handles)
[writefname, writepname] = uiputfile('*','Save SRS data as');
        writepfname = fullfile(writepname, writefname);
        % load data to save in openfile
        fn= getappdata(handles.algsel , 'fn');
        x_pos= getappdata(handles.algsel , 'x_pos');
        x_neg= getappdata(handles.algsel , 'x_neg');
        horlabelt = getappdata(handles.HlabelT , 'horlabelt');
        % WRITE DATA
        fid = fopen(writepfname,'w');
        if horlabelt==1
                fn=1./fn;
                writedata = [fn' x_pos' (abs(x_neg))'];
                fprintf(fid,'  %6s  %12s  %12s\n','Natural Period (s)', 'MaxPositive Response', 'MinPositive Response');    
        else
                writedata = [fn' x_pos' (abs(x_neg))' ];
                fprintf(fid,'  %6s  %12s  %12s\n','Natural Prequency (Hz)', 'MaxPositive Response', 'MinPositive Response');
        end
        fprintf(fid,'  %6g                  %6g                %6g\n',writedata');
        fclose(fid);
guidata(hObject,handles);
% hObject    handle to outputfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in plotPGA.
function plotPGA_Callback(hObject, eventdata, handles)
% load data
fn= getappdata(handles.algsel , 'fn');
x_pos= getappdata(handles.algsel , 'x_pos');
x_neg= getappdata(handles.algsel , 'x_neg');
iunit= getappdata(handles.unitms2 , 'iunit');
horlabelt = getappdata(handles.HlabelT , 'horlabelt');
damp = getappdata(handles.dampsol , 'damp');
Q = getappdata(handles.Qsol , 'Q');
% represent (plot)
if horlabelt==1
    if iunit==0
        handles.ejex=1./fn;
        x_neg=x_neg./9.81; 
        x_pos=x_pos./9.81; 
        plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration (G)');
    else
        handles.ejex=1./fn;
        plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration a_g(m/s^2)');
    end
else
    if iunit==0
      handles.ejex=fn;
      x_neg=x_neg./9.81; 
      x_pos=x_pos./9.81; 
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel ('Peak Acceleration (G)');
    else
      handles.ejex=fn;
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel ('Peak Acceleration a_g(m/s^2)');
    end      
end
title(sprintf('Acceleration Shock Response Spectrum Q=%g',Q));
grid on;
%set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','log');
legend ('positive','negative',2,'Location','best');
guidata(hObject,handles);
% hObject    handle to plotPGA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in HlabelF.
function HlabelF_Callback(hObject, eventdata, handles)
    horlabel = get(handles.HlabelF,'Value');
    fn= getappdata(handles.algsel , 'fn');
    x_pos= getappdata(handles.algsel , 'x_pos');
    x_neg= getappdata(handles.algsel , 'x_neg');
    iunit= getappdata(handles.unitms2 , 'iunit');
    Q = getappdata(handles.Qsol , 'Q');
    if horlabel==1
        horlabelt=0;
        setappdata(handles.HlabelT, 'horlabelt', horlabelt);
    end
    % refresh figure
    if iunit==0
      handles.ejex=fn;
      x_neg=x_neg./9.81; 
      x_pos=x_pos./9.81; 
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel ('Peak Acceleration (G)');
    else
      handles.ejex=fn;
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel ('Peak Acceleration a_g(m/s^2)');
    end  
    title(sprintf('Acceleration Shock Response Spectrum Q=%g',Q));
    grid on;
    legend ('positive','negative',2,'Location','best');
guidata(hObject,handles);
% hObject    handle to HlabelF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of HlabelF
% --- Executes on button press in HlabelT.
function HlabelT_Callback(hObject, eventdata, handles)
    horlabelt = get(handles.HlabelT,'Value');
    fn= getappdata(handles.algsel , 'fn');
    x_pos= getappdata(handles.algsel , 'x_pos');
    x_neg= getappdata(handles.algsel , 'x_neg');
    iunit= getappdata(handles.unitms2 , 'iunit');
    Q = getappdata(handles.Qsol , 'Q');
    setappdata(handles.HlabelT, 'horlabelt', horlabelt);
    % refresh figure
    if iunit==0
        handles.ejex=1./fn;
        x_neg=x_neg./9.81; 
        x_pos=x_pos./9.81; 
        pepe = plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration (G)');
    else
        handles.ejex=1./fn;
        pepe = plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration a_g(m/s^2)');
    end
    title(sprintf('Acceleration Shock Response Spectrum Q=%g',Q));
    grid on;
    legend ('positive','negative',2,'Location','best');
guidata(hObject,handles);
% hObject    handle to HlabelT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of HlabelT
% --- Executes on button press in unitms2.
function unitms2_Callback(hObject, eventdata, handles)
    iunit=get(handles.unitms2,'Value');
    fn= getappdata(handles.algsel , 'fn');
    x_pos= getappdata(handles.algsel , 'x_pos');
    x_neg= getappdata(handles.algsel , 'x_neg');
    Q = getappdata(handles.Qsol , 'Q');
    horlabelt = getappdata(handles.HlabelT , 'horlabelt');
    setappdata(handles.unitms2, 'iunit', iunit);
    % refresh figure
    if horlabelt==1
        handles.ejex=1./fn;
        pepe = plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration a_g(m/s^2)');
    else
      handles.ejex=fn;
      plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
      xlabel('Natural Frequency (Hz)');
      ylabel ('Peak Acceleration a_g(m/s^2)');
    end    
    title(sprintf('Acceleration Shock Response Spectrum Q=%g',Q));
    grid on;
    legend ('positive','negative',2,'Location','best');
guidata(hObject,handles)
% hObject    handle to unitms2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of unitms2
% --- Executes on button press in unitG.
function unitG_Callback(hObject, eventdata, handles)
    aux = get(handles.unitG,'Value');
    fn= getappdata(handles.algsel , 'fn');
    x_pos= getappdata(handles.algsel , 'x_pos');
    x_neg= getappdata(handles.algsel , 'x_neg');
    Q = getappdata(handles.Qsol , 'Q');
    horlabelt = getappdata(handles.HlabelT , 'horlabelt');
    if aux==1
        iunit=1;
        setappdata(handles.unitms2, 'iunit', iunit);
    end
    % refresh figure
    if horlabelt==1
        handles.ejex=1./fn;
        x_neg=x_neg./9.81; 
        x_pos=x_pos./9.81; 
        pepe = plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Period (s)');
        ylabel ('Peak Acceleration (G)');
    else
        handles.ejex=fn;
        x_neg=x_neg./9.81; 
        x_pos=x_pos./9.81; 
        plot(handles.ejex,x_pos,handles.ejex,abs(x_neg),'-.');
        xlabel('Natural Frequency (Hz)');
        ylabel ('Peak Acceleration (G)');
    end
    title(sprintf('Acceleration Shock Response Spectrum Q=%g',Q));
    grid on;
    legend ('positive','negative',2,'Location','best');
guidata(hObject,handles) 
% hObject    handle to unitG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of unitG
%% Aknowledgements
% The author greatly appreciates the code developed by Tom Irvine and
% kindly uploaded as srs.m in 2006. This code is included as one of the
% eligible calculation methods.

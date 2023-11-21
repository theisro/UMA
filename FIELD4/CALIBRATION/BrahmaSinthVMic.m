%*************************************************************************************************
% EM32SinthVMic.m                                                            Lorenzo Chiesi 2009 *
%                                                                                                *
% Questo script esegue la sintesi dei filtri FIR mediante l'inversione della risposta            * 
% all'impulso sferica di un array microfonico permettendo di visualizzare                        *
% ed esportare i risultati                                                                       *
%                                                                                                *
%*************************************************************************************************

if (exist('OCTAVE_VERSION', 'builtin'))
  pkg load signal
  graphics_toolkit('fltk')
  Directivity_2D_Plot = 1;
  Directivity_3D_Plot = 1;
end

addpath('matlab_service_functions')

%Risposta all'impulso sferica utilizzata per la sintesi e e per la
%tracciatura dei diagrammi sferici dei filtri sintetizzati
SIR_File = './BRAHMA_SIR.mat';       %File
SIR_Delay = 0;                  %Numero di campioni da saltare
SIR_N = 1024;                   %Numero di campioni da utilizzare


%Risposta all'impulso planare orizzontale e verticale utilizzata per tracciare i
%diagrammi polari dei filtri sintetizzati
HIR_File = './BRAHMA_HIR.mat';       %File
HIR_Delay = 0;                  %Numero di campioni da saltare
HIR_N = 1024;                   %Numero di campioni da utilizzare
VIR_File = '';       %File
VIR_Delay = 0;                  %Numero di campioni da saltare
VIR_N = 1024;                   %Numero di campioni da utilizzare

%Parametri di regolarizzazione
fl = 40;                %Frequenza di taglio inferiore [Hz]
fh = 16000;             %Frequenza di taglio superiore [Hz]
bt = 0.3;               %Banda di transizione [ottave]
Bol = 1;                %Parametro di regolarizzazione fuori-banda h
Bi = 0.001;              %Parametro di regolarizzazione in-banda
Boh = 1;                %Parametro di regolarizzazione fuori-banda l

%Lunghezza dei filtri sintetizzati (FIR_N >= SIR_N)
FIR_N = 2048;

%Abilitazione grafici
SIR_Plot                = 0;
HIR_Plot                = 0;
VIR_Plot                = 0;
FIR_Plot                = 0;
if (exist('OCTAVE_VERSION', 'builtin') == 0)
    Directivity_3D_Plot     = 1;
    Directivity_2D_Plot     = 1;
end

%Esportazione del termine indipendente dalla direttivitŕ
R_File = ''; %'R_matrix_1024.dat';


%Esportazione dei filtri FIR sintetizzati
FIR_File = 'BRAHMA_A2BFORMAT.wav';         %File WAV in cui esportare i filtri FIR sintetizzati
BFIR_Directory = '';%'FIR';         %Cartella in cui esportare i filtri sintetizzati in formato "brutefir"




%*******************************************************************************************
% Caricamento delle risposta all'impulso sferica "sir" (Dir,Mic,SIR_N) 
%*******************************************************************************************


if not(exist('sir','var'))
    
    %Caricamento delle risposta all'impulso sferica
    load(SIR_File);
    Dir   = size(sir,1);
    Mic   = size(sir,2);
    Fs    = sir_Fs;
    
    %Ritaglio delle risposta all'impulso sferica
    sir = sir(:,:,SIR_Delay+1:SIR_Delay+SIR_N);
    
    if (SIR_Plot)
        %Grafico risposta all'impulso sferica        
        figure('Name','Risposta all''impulso sferica (verificare ritaglio)','NumberTitle','off')
        plot(squeeze(sir(1,:,:))')  
        axis([1 SIR_N -1 +1])        
    end

end



%*******************************************************************************************
% Calcolo della direttivitŕ dei microfoni virtuali "A" (Position, VirtualMic)
%*******************************************************************************************

%Numero di microfoni virtuali
VirtualMic = 4;     

%Rotazione globale del sistema di riferimento
pan = 0;
tilt = 0;


%Calcola il valore della target function per ogni direzione
A = zeros(Dir,VirtualMic);
for Dir_i = 1:Dir

    %Coordinate sferiche della direzione corrente
    theta_dir = sir_angle(Dir_i,1);
    phi_dir   = sir_angle(Dir_i,2);
    
    %Rotazione del sistema di riferimento
    [theta,phi] = sph_rotation(theta_dir,phi_dir,pan,tilt);
        
     
    
%     %********** Pressure/Velocity Microphone ****************
%     
%     %Pressure/Velocity Mix
%     %P = 1;         G = 0;          %Omnidirectional
%     %P = 0.75;      G = 0.25;       %Subcardioid
%     P = 0.5;       G = 0.5;        %Cardioid
%     %P = 0.333;     G = 0.666;      %Supercardioid
%     %P = 0.25;      G = 0.75;       %Hypercardioid
%     %P = 0;         G = 1;          %Bidirectional
%     
%     order = 1;                      %Microphone order
% 
%     fod = P + G .* cos(theta) .* cos(phi);        %First order directicity
%     A(Dir_i,1) = sign(fod).*abs(fod).^order;      %Higher order directivity

%     
%      %******************* Anello orizzontale di microfoni ************************* 
%      order = 2;
%      P = 0.5;       
%      G = 0.5;
%      VMicDir = [ (0:2*pi/VirtualMic:2*pi*((VirtualMic-1)/VirtualMic))' , zeros(VirtualMic,1) ];
%      for i = 1:VirtualMic
%          [phi_t,theta_t] = sph_rotation(theta, phi, VMicDir(i,1), VMicDir(i,2));    
%          A(Dir_i,i) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%      end

    
%     %******************* Stereo XY *****************************
%        
%     A(Dir_i,1) = 0.5 + 0.5 * cos(theta-pi/4) * cos(phi);      %L
%     A(Dir_i,2) = 0.5 + 0.5 * cos(theta+pi/4) * cos(phi);      %R       



%     %******************* Surround 5.1 *************************
%        
%     A(Dir_i,1) = 0.33 + 0.66 * cos(theta-pi/3) * cos(phi);      %L
%     A(Dir_i,2) = 0.33 + 0.66 * cos(theta+pi/3) * cos(phi);      %R
%     A(Dir_i,3) = 0.25 + 0.75 * cos(theta) * cos(phi);           %C
%     A(Dir_i,4) = 0.5  + 0.5  * cos(theta-3/4*pi) * cos(phi);    %BL
%     A(Dir_i,5) = 0.5  + 0.5  * cos(theta+3/4*pi) * cos(phi);    %BR
%     A(Dir_i,6) = 1;                                             %S



%    ************************** Ambisonic Armonics *************************
     
    A(Dir_i,1) = sqrt(1/2);                                             %W
    A(Dir_i,2) = cos(theta)*cos(phi);                                   %X
    A(Dir_i,3) = sin(theta)*cos(phi);                                   %Y
    A(Dir_i,4) = sin(phi);                                              %Z 
          
%     A(Dir_i,5) = 1.5*sin(phi)^2 - 0.5;                                  %R
%     A(Dir_i,6) = cos(theta)*sin(2*phi);                                 %S
%     A(Dir_i,7) = sin(theta)*sin(2*phi);                                 %T
%     A(Dir_i,8) = cos(2*theta)*cos(phi)^2;                               %U
%     A(Dir_i,9) = sin(2*theta)*cos(phi)^2;                               %V
%        
%     A(Dir_i,10) = sin(phi)*(5*sin(phi)^2 - 3)/2;                        %K
%     A(Dir_i,11) = sqrt(135/256)*cos(theta)*cos(phi)*(5*sin(phi)^2 - 1); %L
%     A(Dir_i,12) = sqrt(135/256)*sin(theta)*cos(phi)*(5*sin(phi)^2 - 1); %M
%     A(Dir_i,13) = sqrt(27/4)*cos(2*theta)*sin(phi)*cos(phi)^2;          %N
%     A(Dir_i,14) = sqrt(27/4)*sin(2*theta)*sin(phi)*cos(phi)^2;          %O
%     A(Dir_i,15) = cos(3*theta)*cos(phi)^3;                              %P
%     A(Dir_i,16) = sin(3*theta)*cos(phi)^3;                              %Q
%  


%     %********************* Sala CDM ***************************
% 
%     %P = 0.75;      G = 0.25;       %Subcardioid
%     P = 0.5;       G = 0.5;        %Cardioid
%     %P = 0.333;     G = 0.666;      %Supercardioid
%     %P = 0.25;      G = 0.75;       %Hypercardioid
%         
%     order = 3;                      %Microphone order
%         
%     [theta_t,phi_t] = sph_rotation( theta, phi,   0/180*pi,     0 );    A(Dir_i,1) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi,  45/180*pi,     0 );    A(Dir_i,2) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi,  90/180*pi,     0 );    A(Dir_i,3) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 135/180*pi,     0 );    A(Dir_i,4) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 180/180*pi,     0 );    A(Dir_i,5) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 225/180*pi,     0 );    A(Dir_i,6) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 270/180*pi,     0 );    A(Dir_i,7) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 315/180*pi,     0 );    A(Dir_i,8) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi,   0/180*pi,  pi/4 );    A(Dir_i,9) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi,  90/180*pi,  pi/4 );    A(Dir_i,10) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 180/180*pi,  pi/4 );    A(Dir_i,11) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 270/180*pi,  pi/4 );    A(Dir_i,12) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi,   0/180*pi, -pi/4 );    A(Dir_i,13) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi,  90/180*pi, -pi/4 );    A(Dir_i,14) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 180/180*pi, -pi/4 );    A(Dir_i,15) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%     [theta_t,phi_t] = sph_rotation( theta, phi, 270/180*pi, -pi/4 );    A(Dir_i,16) = hod(P + G * cos(theta_t) * cos(phi_t) ,order); 
%      


%      %******************* Rumore Farina (24 Mic) *************************
%      order = 5;
%      P = 0.5;       G = 0.5;        %Cardioid
%      for i = 1:24
%          fod = P + G .* cos(theta+((i-1)*(2*pi)/24)) .* cos(phi);        %First order directicity
%          A(Dir_i,i) = sign(fod).*abs(fod).^order;      %Higher order directivity
%      end
 
%     %************************** Omni + 8 Figure *************************
%      
%      A(Dir_i,1) = 1;                       %O
%      A(Dir_i,2) = sin(theta)*cos(phi);     %8
% 

%      %******************* La Scala Farina (32 Mic) *************************
%      GenerateEM32Dir;
%      order = 4;
%      P = 0.5;       
%      G = 0.5;
%      for i = 1:VirtualMic
%          [phi_t,theta_t] = sph_rotation(theta, phi, EM32Dir(i,2), EM32Dir(i,3));    
%          A(Dir_i,i) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%      end


%      %******************* 122 Mic *************************
%      VMicDir = GenerateUniformTestDir(4);
%      order = 16;
%      P = 0.5;       
%      G = 0.5;
%      for i = 1:VirtualMic
%          [phi_t,theta_t] = sph_rotation(theta, phi, VMicDir(i,1), VMicDir(i,2));    
%          A(Dir_i,i) = hod(P + G * cos(theta_t) * cos(phi_t) ,order);
%      end

end




%*******************************************************************************************
%                           Inversione delle risposte all'impulso
%*******************************************************************************************

W = ones(Dir,1,'single');
[h,R] = invert_kirkeby(sir,A,W,Fs,fl,fh,bt,Bol,Boh,Bi,FIR_N);

%r = ifft( cat(3,R, flipdim( conj(R(:,:,2:end-1)) ,3)) ,[],3);
%figure
%spherical_directivity2(Fs,sir,sir_angle,r(:,167,:),sir_angle(167,1),sir_angle(167,2))

%*******************************************************************************************
%                                Visualizzazione dei risultati
%*******************************************************************************************

%Visualizza coefficienti dei filtri FIR (Limitato al primo microfono virtuale) 
if (FIR_Plot)
    figure('Name','Filtri FIR sintetizzati (VMic 1) (verificare la convergenza)','NumberTitle','off')
    plot(squeeze(h(:,1,:))')
end




%Plotta la direttivitŕ 3D
%if (Directivity_3D_Plot)
  %  spherical_directivity2(Fs,sir,sir_angle,h,0,0)
%end




%Plotta la direttivitŕ 2D
if (Directivity_2D_Plot)
   
    %Carica e ritaglia la risposta all'impulso orizzontale
    if not(exist('hir','var'))
        
        load(HIR_File);
        hir = hir(:,:,HIR_Delay+1:HIR_Delay+HIR_N);
        
        if (HIR_Plot)     
            %Grafico risposta all'impulso sferica        
            figure('Name','Risposta all''impulso orizzontale (verificare ritaglio)','NumberTitle','off')
            plot(squeeze(hir(1,:,:))')  
            axis([1 HIR_N -1 +1])        
        end
    end
    
    if (Directivity_2D_Plot)
        plane_directivity(Fs,hir,h,0) 
    end
    
%     %Carica e ritaglia la risposta all'impulso verticale
%     if not(exist('vir','var'))
%         
%         load(VIR_File);
%         vir = vir(:,:,HIR_Delay+1:HIR_Delay+HIR_N);
%         
%         if (VIR_Plot)     
%             %Grafico risposta all'impulso sferica        
%             figure('Name','Risposta all''impulso verticale (verificare ritaglio)','NumberTitle','off')
%             plot(squeeze(vir(1,:,:))')  
%             axis([1 VIR_N -1 +1])        
%         end
%     end
%     
%     plane_directivity(Fs,[vir(10:36,:,:);vir(1:9,:,:)],h)
        
end

%*******************************************************************************************
%                                PATCH PER M+S ZOOM H2N on channel 3 and 4 
%*******************************************************************************************

hms = h;
h(3,:,:) = hms(3,:,:) + hms(4,:,:);
h(4,:,:) = hms(3,:,:) - hms(4,:,:);




%*******************************************************************************************
%                                Esportazione dei risultati
%*******************************************************************************************

%Memorizza il termine indipendente dalla direttivitŕ
%per l'uso con ilsistema di elaborazione real-time
if (R_File)
    fprintf('Esportazione termine indiendente dalla direttivitŕ: %s\n',R_File)
    store_matrix(sprintf('%s.dat',R_File),R);
    R_angle = sir_angle(:,1)+1i*sir_angle(:,2);
    store_matrix(sprintf('%s.dir',R_File),R_angle);
    %clear R
end
    
       
    
%Memorizza il filtro inverso sintetizzato
if (FIR_File)
   fprintf('Esportazione filtri FIR compatibili con XVolver: %s\n',FIR_File)
   store_fir(FIR_File, h, Fs)
end

%Memorizza filtro per brutefir
if (BFIR_Directory)
    fprintf('Esportazione filtri FIR compatibili con BruteFIR nella cartella: %s\n',BFIR_Directory)
    store_brutefir(BFIR_Directory,h)
end


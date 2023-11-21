%*************************************************************************************************
%CylPackHorizontalIR.m                                                      Lorenzo Chiesi 2013 *
%                                                                                               *
%                                                                                               *
%*************************************************************************************************

close all

Directivity_3D_Plot = 1;
Directivity_2D_Plot = 1;
if (exist('OCTAVE_VERSION', 'builtin'))
  pkg load signal
  graphics_toolkit('fltk')
  Directivity_2D_Plot = 1;
  Directivity_3D_Plot = 1;
end

addpath('matlab_service_functions')

%Parametri impostati per array EM32

Mic           = 4;           %Numero di microfoni
Dir = 8;                     %Numero di direzioni misurate
StepAngle = 360/Dir;
sound_speed   = 343;         %Velocitŕ del suono [m/s]

N             = 1024;        %Lunghezza delle IR   

rhir_angle = [0 45 90 135 180 225 270 315]/180*pi;

%*************************************************************************************************
%                         CARICAMENTO E EQUALIZZAZIONE DELLE IR
%*************************************************************************************************


%Parametri di regolarizzazione
fl = 40;                %Frequenza di taglio inferiore [Hz]
fh = 16000;             %Frequenza di taglio superiore [Hz]
bt = 0.3;               %Banda di transizione [ottave]
Bol = 1;                %Parametro di regolarizzazione fuori-banda h
Bi = 0.001;             %Parametro di regolarizzazione in-banda
Boh = 1;                %Parametro di regolarizzazione fuori-banda l

%Load and win calibration IR
[calIR,Fs] = audioread('./calibrate.wav');
calIR = TrimWinData( calIR(:,1)', 0, N/2, 'hann', 240, 1, 64, 'calibration_trim_plot' );
%return;

%Compute inverse filter
postfilter = invert_ir(calIR,Fs,fl,fh,bt,Bol,Boh,Bi,N);

figure
irplot(postfilter)

%definisci matrice delle IR (AzimuthDir x Mic x N)
ir = zeros(Dir,Mic,N,'single');

for d = 1:Dir
    
    %Carica WAV stereo relativi alla direzione selezionata
    [FrontIR,~] = audioread(sprintf('./front/%d.wav',(d-1)*360/Dir));
    [RearIR,Fs] = audioread(sprintf('./rear/%d.wav',(d-1)*360/Dir));
    
    %Store to matrix
    ir(d,1:2,:) = FrontIR(1:N,:)';
    ir(d,3:4,:) = RearIR(1:N,:)';
end

%Apply equalization
ireq = zeros(Dir,Mic,N,'single');
for  d=1:Dir
    for m=1:Mic

        res = fd_conv(squeeze(ir(d,m,:)),postfilter);
        ireq(d,m,:) = res(1:N);
            
    end
end

hir = TrimWinData(ireq, 300, N, 'hann', 400, 64, 64, 'ir_trim_plot');
%return;


%*************************************************************************************************
%         ANALISI DELLA DISTANZA TRA L'ARRAY E LA SORGENTE PER OGNI DIREZIONE DI MISURA
%*************************************************************************************************

ft = fd_flyingtime_estimation(hir,250/Fs,8000/Fs);
%ft = flyingtime_estimation2(hir);

skew = mean(ft,2);
skew = skew - mean(skew);


%Utilizzando un approccio numerico ai minimi quadrati interpoliamo la variazione della distanza
%con una funzione sinusoidale del tipo: a * cos(table_angle + b) + c
norm(1) = (max(skew)-min(skew))/2;
norm(2) = 2*pi;
norm(3) = mean(skew);
f = @(var) sum( ( var(1)*norm(1) * cos(rhir_angle + var(2)*norm(2)) + var(3)*norm(3) - skew' ).^2 );
options = optimset('TolFun',1e-6);
var = fminsearch(f, [1,0,1], options);
interp_skew = ( var(1)*norm(1) * cos(rhir_angle + var(2)*norm(2)) + var(3)*norm(3) );

show_angle = linspace(0,2*pi,1000);
show_interp_skew = ( var(1)*norm(1) * cos(show_angle + var(2)*norm(2)) + var(3)*norm(3) );

fprintf('\n\nOscillazione asse "theta" :\n')
fprintf('Ampiezza: %f [sample], %f [mm]\n',var(1)*norm(1),var(1)*norm(1)/Fs * sound_speed * 1e3);
fprintf('Fase: %f [rad], %f [deg]\n\n',var(2)*norm(2),var(2)*norm(2)/pi*180);


%Converte l'errore di posizione in mm e lo visualizza in un grafico al variare dell'angolo della tavola rotante.
figure
table_angle = 0:StepAngle:StepAngle*(Dir-1);
plot(table_angle,(skew/Fs * sound_speed * 1e3)','*',show_angle/pi*180,(show_interp_skew/Fs * sound_speed * 1e3),'-')
title('Spostamento del centro acustico dell''array microfonico durante la rotazione');
ylabel('Spostamento dalla posizione intermedia [mm]');
xlabel('Posizione della tavola rotante [°]');


%Compensazione dell'errore sinusoidale
unskew = - interp_skew;                

%Compensazione totale dell'errore
%unskew = - skew;       

%Nessuna compensazione
%unskew = zeros(size(skew,1),size(skew,2));


%*************************************************************************************************
%             APPLCAZIONE DELLE TRASLAZIONI TEMPORALI NEL DOMINIO DELLA FREQUENZA
%*************************************************************************************************
disp('Applicazione delle traslazioni temporali nel dominio della frequenza...')

%Applica le traslazioni temporali ad ogni IR
for Azimuth_i = 1:Dir
    for Mic_i = 1:Mic
        hir(Azimuth_i,Mic_i,:) = fd_time_shift(squeeze(hir(Azimuth_i,Mic_i,:))',unskew(Azimuth_i));
    end
end




%*************************************************************************************************
%                         HIR2SIR
%*************************************************************************************************


%Extract 8 ir from hir
rhir = hir;
rhir_dir = 8;
rhir_angle = [0 45 90 135 180 225 270 315]/180*pi;


FIR_N = 1024;


FitPoint = [
 1 1 1 1 1 0 1 1 
 1 1 1 0 1 1 1 1 
 1 1 1 1 1 1 1 0 
 1 0 1 1 1 1 1 1 
]';


%Define Microphone capsules directions
% http://pcfarina.eng.unipr.it/Public/B-format/A2B-conversion/A2B.htm
%L/R (Left,Right); F/B (Front,Back); D/U (Down,Up)
CartMicDir = [     
    1,  1,  1;  %1 LFU
    1, -1, -1;  %2 RFD
   -1,  1, -1;  %3 LBD
   -1, -1,  1;  %4 RBU
];
%Convert to spherical
[theta,phi,r] = cart2sph(CartMicDir(:,1),CartMicDir(:,2),CartMicDir(:,3));
MicDir = [theta,phi];


%Plot IR loaded
 
%Evaluate flingtime to extract array radius
%ft = flyingtime_estimation2(double(rhir));
ft = fd_flyingtime_estimation(double(rhir),250/Fs,8000/Fs);
ftm = round(mean(mean(ft)));
ft = ft - mean(mean(ft));

mic_r = zeros(3,1);
mic_t = zeros(3,1);
mic_o = zeros(3,1);
figure
for m=1:4

    MaskSel = logical(FitPoint(:,m));
    
    ft2D = -ft(MaskSel,m);
    rhir_angle2 = rhir_angle(MaskSel);
    
    
    %Utilizzando un approccio numerico ai minimi quadrati interpoliamo la variazione della distanza
    %con una funzione sinusoidale del tipo: a * cos(theta + b) + c
    f = @(var) sum( ( var(1) * cos(rhir_angle2 - MicDir(m,1) + var(2)) + var(3) - ft2D' ).^2 );
    var = fminsearch(f, [1,0,0]);
    fit_angle = linspace(0,2*pi,1000);
    fit_ft = var(1)*cos(fit_angle-MicDir(m,1)+var(2))+var(3);
    
    if var(1)<0
        var(1) = -var(1);
        var(2) = var(2)+pi;
        var(2) = atan2(sin(var(2)),cos(var(2)));
    end

    subplot(4,2,(m-1)*2+1)
    irplot(squeeze(rhir(:,m,ftm-8:ftm+8)))
    subplot(4,2,(m-1)*2+1+1)
    plot(rhir_angle2/pi*180,ft2D,'r*',fit_angle/pi*180,fit_ft,'b')

    mic_r(m) = var(1);
    mic_t(m) = var(2);
    mic_o(m) = var(3);
    
end

dpa_r = mean(mic_r./cos(MicDir(:,2)));
dpa_t = mean(mic_t);
 
%Applica le traslazioni temporali ad ogni IR
for d = 1:rhir_dir
    for m = 1:4
        rhir(d,m,:) =  fd_time_shift( rhir(d,m,:) , dpa_r*cos(rhir_angle(d)-MicDir(m,1)+dpa_t) );
    end
end




%From each Mic, measured in 8 different directions I want to compute:
%1) Pressure component, 2) Pressure Gradient component
PGir = zeros(4,2,FIR_N);
figure
for m = 1:1:4
 
    %Compute angle between Mic axis direction and each measurement directions
    ang_dist = sph_dist( MicDir(m,1),MicDir(m,2) , rhir_angle,0 );
    

    %1) 1*Pir + cos(a1)*Gir = ir1
    %2) 1*Pir + cos(a2)*Gir = ir2
    %   ...
    %D) 1*Pir + cos(aD)*Gir = irD
    
    %[Dx2] [2x1] [Dx1]
    %  C  *  H  =  A
    
     C_matrix = [ ones(rhir_dir,1) , cos(ang_dist)'];
     A_matrix = fft( rhir(:,m,:) , FIR_N, 3);
     H_matrix = zeros(2,1,FIR_N); 

     for k = 2 : (FIR_N/2)
        H_matrix(:,:,k) = ( (C_matrix(:,:)')*C_matrix(:,:) ) \ ( C_matrix(:,:)' * A_matrix(:,:,k)  );   
     end
     
    for i=1:2     
        H_matrix(i,:,FIR_N/2+1:FIR_N) = [0; conj( flipud( squeeze(H_matrix(i,:,2:FIR_N/2)) ) )];
    end

    PGir(m,:,:) = permute( ifft(H_matrix,FIR_N,3) , [2 1 3] ); 
    
    subplot(4,1,m)
    irplot(squeeze(PGir(m,:,ftm-16:ftm+64)))
end







%****************************************
%Create virtual hir matrix
%****************************************
hir_dir = 72;
hir_angle = 0:2*pi/hir_dir:2*pi*(1-1/hir_dir); 
hir = zeros(hir_dir,4,FIR_N);
for d=1:1:hir_dir
    for m=1:1:4
        ang_dist = sph_dist( hir_angle(d),0 , MicDir(m,1),MicDir(m,2) );
        hir(d,m,:)= PGir(m,1,:) + PGir(m,2,:)*cos(ang_dist);
        hir(d,m,:) = fd_time_shift( hir(d,m,:) , -dpa_r*cos(ang_dist) );
    end
end

if (Directivity_2D_Plot)
    plane_directivity(Fs,hir)
end

%Export virtual sir matrix
hir_Fs = Fs;
save -v7 'BRAHMA_HIR.mat' hir hir_angle hir_Fs



%****************************************
%Create virtual sir matrix
%****************************************


%Generate uniform sir test direction pattern
sir_angle = GenerateUniformTestDir(4);
Dir = length(sir_angle);

%Create sir matrix [Dir,Mic,N]
sir = zeros(Dir,4,FIR_N);
for d=1:1:Dir
    for m=1:1:4
       ang_dist = sph_dist( sir_angle(d,1),sir_angle(d,2) , MicDir(m,1),MicDir(m,2) );
       sir(d,m,:) = PGir(m,1,:) + PGir(m,2,:)*cos(ang_dist); 
       sir(d,m,:) = fd_time_shift( sir(d,m,:) , -dpa_r*cos(ang_dist) );
    end
end

%Plot resulting 3D pattern
%if (Directivity_3D_Plot)
 %   spherical_directivity2(Fs,sir,sir_angle)
%end

%Export virtual sir matrix
sir_Fs = Fs;
save -v7 'BRAHMA_SIR.mat' sir sir_angle sir_Fs


%%********************************************************************************
%*               STORE FIR MATRIX IN A MULTICHANNEL WAV FILE                     *
%*          by Lorenzo Chiesi 2009 - lorenzo.chiesi(at)gmail.com                 *
%*********************************************************************************
%
% store_fir(fir_file,fir_matrix,fs,gain)
%
% Input:  fir_file   : WAV file to store the FIR matrix
%         fir_matrix : Matrix of FIR filter (In x Out x N)
%         fs         : Sampling frequency
%         dB_gain    : dB gain. If omitted normalization is performed


function store_fir(fir_file,fir_matrix,fs,dB_gain)

%Calcola parametri
in = size(fir_matrix,1);    %Numero di ingressi della matrice di filtri
out = size(fir_matrix,2);   %Numero di uscite della matrice di filtri
N = size(fir_matrix,3);     %Lunghezza dei filtri


fir_max = squeeze(max(max(max(abs(fir_matrix)))));     %Calcola il valore massimo

if not( exist('dB_gain','var') )
    gain = 1/fir_max;                                  %Calcola il coefficiente del guadagno
    dB_gain = 20*log10(gain);
    fprintf('Store FIR filter in "%s". (Normalization gain: %0.3f dB)\n',fir_file, dB_gain);
else
    gain = 10^(dB_gain/20);
    fprintf('Store FIR filter in "%s". (Imposed gain: %0.3f dB)\n',fir_file, dB_gain);
    if (abs(fir_max*gain) > 1)
        fprintf('!!WARNING!! WAV FILE CLIPPED!! Max: %f\n',fir_max*gain);
    end
end

%Prealloca memoria matrice per output su file (float 32 bit)
wav = zeros(in*N,out,'single');

%Riarrangia matrice del filtro nel formato adatto alla scrittura su file
for in_i = 1:in
    wav((in_i-1)*N+1:in_i*N,:) = (squeeze(fir_matrix(in_i,:,:))').*gain;
end

clear fir_matrix

%Scrivi su file
%wavwrite(wav,fs,32,fir_file);
%This add metedata to support X-Volver
audiowrite(fir_file,wav,fs,'BitsPerSample',32,'Comment','XvMulti');

end


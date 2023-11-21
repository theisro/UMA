%*************************************************************************************************
% invert_kirkeby.m                                                           Lorenzo Chiesi 2009 *
%                                                                                                *
% Questa funzione implementa il processo di sintesi dei filtri FIR necessari per l'elaborazione  * 
% del segnale di un array microfonico mediante l'inversione nel dominio della frequenza          *  
% delle risposte all'impulso misurate sperimentalmente                                           *                      
%                                                                                                *
%*************************************************************************************************
%
% [h_matrix,R_matrix] = invert_kirkeby(c_matrix,A_matrix,w_vector,Fs,fl,fh,bt,Bol,Bi,Boh,Nh);
%
% c_matrix  = Matrice delle risposte all'impulso misurate (Dir x Mic x Nc)
% A_matrix  = Matrice dei guadagni desiderati (Dir x VirtualMic )
% w_vector  = Vettore dei coefficienti di peso delle posizioni di misura (Dir)
% fl        = Frequenza di taglio inferiore [Hz]
% fh        = Frequenza di taglio superiore [Hz]
% bt        = Banda di transizione [ottave]
% Bol       = Parametro di regolarizzazione fuori-banda in BF
% Bi        = Parametro di regolarizzazione in-banda
% Boh       = Parametro di regolarizzazione fuori-banda in HF 
% Nh         = Lunghezza dei filtri FIR sintetizzati (Nh > Nc)
%
% h_matrix  = Matrice dei filtri inversi sintetizzati (Mic x VirtualMic x Nh)
% R_matrix  = Termine indipendente dalla direttivitŕ (Mic x Dir x Nh/2+1)

function [h_matrix,R_matrix] = invert_kirkeby(c_matrix,A_matrix,W_vector,Fs,fl,fh,bt,Bol,Boh,Bi,Nh)

%display(sprintf('Inversione delle risposte all''impulso secondo Kirkeby:'))
    
Dir = size(c_matrix,1);
Mic = size(c_matrix,2);
Nc = size(c_matrix,3);
VirtualMic = size(A_matrix,2);




%********************* Calcolo di C(k) ********************

%display(sprintf('\t- Trasformazione nel dominio della frequenza...'))

%Trasformata delle risposte all'impulso sul numero di campioni desideratiper i filtri
C_matrix = fft(c_matrix,Nh,3);
%Mantieni solo le prime Nh/2+1 righe spettrali per risparmiare memoria 
C_matrix = C_matrix(:,:,1:Nh/2+1);
clear c_matrix




%********************* Calcolo di B(k) *********************

%display(sprintf('\t- Calcolo del parametro di regolarizzazione...'))

% Profilo del parametro Beta al variare della frequenza
%     A
%  Boh|                                   -------
%  Bol+------                            / 
%     |       \                         /
%     |         \                     /
%     |           \                 /
%  Bi +             ---------------
%     +-----+---+---+-------------+---+---+------+---->
%     0    fl1  fl fl2           fh1  fh fh2    Fs/2
                
%Calcolo delle frequenze di transizione
fl1 = fl / (2^(bt/2));
fl2 = fl * (2^(bt/2));
fh1 = fh / (2^(bt/2));
fh2 = fh * (2^(bt/2));

%Quantizzazione dele frequenze di transizione
df = Fs/Nh;
kl1 = round(fl1/df);
kl2 = round(fl2/df)+1;
kh1 = round(fh1/df);
kh2 = min(round(fh2/df)+1,Nh/2);

%Creazione del profilo di B in funzione della frequenza
%B č quindi un vettore con Nh/2+1 elementi (basta rappresentarlo per le frequenze positive)
B = [ones(1,kl1).*Bol, linspace(Bol,Bi,kl2-kl1), ones(1,kh1-kl2).*Bi, linspace(Bi,Boh,kh2-kh1), ones(1,(Nh/2+1)-kh2).*Boh ];





%********************* Applica coefficienti di peso *********************

%Verifica se i coefficienti sono diversi da 1
if (mean(W_vector == ones(size(W_vector,1),size(W_vector,2))) ~= 1) 
    
    %display(sprintf('\t- Applicazione coefficienti di peso...'))
    
    %Applica coefficienti di peso
    for Dir_i = 1:Dir 
        C_matrix(Dir_i,:,:) = sqrt(W_vector(Dir_i)).*C_matrix(Dir_i,:,:);
        A_matrix(Dir_i) = sqrt(W_vector(Dir_i)).*A_matrix(Dir_i);
    end
    
end





%********************* Inversione nel dominio della frequenza *********************

%NOTA: I segnali audio nel dominio del tempo sono ovviamente reali, i loro spettri godono quindi di simmetria hermitiana.
%      Possiamo quindi eseguire il calcolo del filtro inverso H solo per le prime N/2+1 righe spettrali dell'fft 
%      (che in Matlab corrispondono alle frequenze positive) e ottenere le restanti per simmetria.
% in realtŕ la frequenza 1 rappresenta la DC e la frequenza Nh/2+1 rappresenta la Fs e non possono essere rifasate...quindi vengono fissate a 0

%display(sprintf('\t- Inversione...'))

%Inizializza matrici
H_matrix = zeros(Mic,VirtualMic,Nh,'single'); 
R_matrix = zeros(Mic,Dir,Nh/2+1,'single');
%warning('OFF','MATLAB:nearlySingularMatrix')

S = Nh/2;    %Ritardo di causalizzazione

%Calcola R e H per ogni riga spettrale
%(Salta righe spettrali DC e Fs impostate a 0 durante l'inizializzazione)
for k = 2 : (Nh/2)
    
%     %Forma ols calcolata da matlab (senza regolarizzazione) 
%     H_matrix(:,:,k) = C_matrix(:,:,k) \ ((1/2)*A_matrix.*exp(-j*2*pi*(k-1)*((S+1)-1)/Nh));  

%     %Forma con regolarizzazione
%     H_matrix(:,:,k) = ( (C_matrix(:,:,k)')*C_matrix(:,:,k) + B(k).*eye(Mic)) \ ( (C_matrix(:,:,k)') * ((1/2)*A_matrix.*exp(-j*2*pi*(k-1)*((S+1)-1)/Nh)));  
   
%     %Forma adatta all'esportazione
%     R_matrix(:,:,k) = inv( (C_matrix(:,:,k)')*C_matrix(:,:,k) + B(k).*eye(Mic)) * (C_matrix(:,:,k)').*(1/2)*exp(-j*2*pi*(k-1)*((S+1)-1)/Nh);
%     H_matrix(:,:,k) = ( R_matrix(:,:,k) ) * (A_matrix);  
    
    %Forma adatta all'esportazione (Avoid use of inv that is deprecated by MATLAB)
    R_matrix(:,:,k) = ( (C_matrix(:,:,k)')*C_matrix(:,:,k) + B(k).*eye(Mic) ) \ (C_matrix(:,:,k)') .* (1/2)*exp(-j*2*pi*(k-1)*((S+1)-1)/Nh);
    H_matrix(:,:,k) = ( R_matrix(:,:,k) ) * (A_matrix);  
    

end

%Elimina matrici non piů utilizzate per risparmiare memoria
clear C_matrix

%Rendi hermitiana la matrice H 
for Mic_i=1:Mic
    for VirtualMic_i=1:VirtualMic       
        H_matrix(Mic_i,VirtualMic_i,Nh/2+1:Nh) = [0; conj( flipud( squeeze(H_matrix(Mic_i,VirtualMic_i,2:Nh/2)) ) )];
    end
end





%********************* Calcolo di h(n) *********************


%display(sprintf('\t- Antitrasformazione delle funzioni di trasferimento...'))
h_matrix = ifft(H_matrix,Nh,3);        
clear H_matrix


%display(sprintf('\t- Finestratura dei filtri FIR...'))

%Costruisci la finestra
twin = triang(Nh/2);
window = [twin(1:floor(Nh/4));ones(floor(Nh/2),1);twin(floor(Nh/4)+1:floor(Nh/2))];
window = [window ; zeros(size(h_matrix,3)-size(window,1),1)];

%Applica la finestra
for Mic_i=1:Mic
    for VirtualMic_i=1:VirtualMic
        h_matrix(Mic_i,VirtualMic_i,:) = squeeze(h_matrix(Mic_i,VirtualMic_i,:)).*window;
    end
end



end




%*************************************************************************************************
% invert_ir.m                                                                Lorenzo Chiesi 2011 *
%                                                                                                *
% Single channel kyrkeby inversion                                                               *                      
%                                                                                                *
%*************************************************************************************************
%
% fir_matrix = invert_ir(ir_matrix,Fs,fl,fh,bt,Bol,Boh,Bi,No);
%
% ir_matrix   = Matrice delle risposte all'impulso misurate (M x Ni)
% fl          = Frequenza di taglio inferiore [Hz]
% fh          = Frequenza di taglio superiore [Hz]
% bt          = Banda di transizione [ottave]
% Bol         = Parametro di regolarizzazione fuori-banda in BF
% Boh         = Parametro di regolarizzazione fuori-banda in HF
% Bi          = Parametro di regolarizzazione in-banda
% No          = Lunghezza dei filtri FIR sintetizzati (No > Ni)
%
% fir_matrix  = Matrice dei filtri inversi sintetizzati (M x No)

function fir_matrix = invert_ir(ir_matrix,Fs,fl,fh,bt,Bol,Boh,Bi,No)
    
M = size(ir_matrix,1);

%Trasformata delle risposte all'impulso sul numero di campioni desiderati per i filtri
IR_matrix = fft(ir_matrix,No,2);
%Mantieni solo le prime Nh/2+1 righe spettrali per risparmiare memoria 
IR_matrix = IR_matrix(:,1:No/2+1);
clear ir_matrix


%********************* Calcolo di B(k) *********************

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
df = Fs/No;
kl1 = round(fl1/df);
kl2 = round(fl2/df)+1;
kh1 = round(fh1/df);
kh2 = min(round(fh2/df)+1,No/2);

%Creazione del profilo di B in funzione della frequenza
%B č quindi un vettore con No/2+1 elementi (basta rappresentarlo per le frequenze positive)
B = [ones(1,kl1).*Bol, linspace(Bol,Bi,kl2-kl1), ones(1,kh1-kl2).*Bi, linspace(Bi,Boh,kh2-kh1), ones(1,(No/2+1)-kh2).*Boh ];


%********************* Inversione nel dominio della frequenza *********************

%NOTA: I segnali audio nel dominio del tempo sono ovviamente reali, i loro spettri godono quindi di simmetria hermitiana.
%      Possiamo quindi eseguire il calcolo del filtro inverso H solo per le prime N/2+1 righe spettrali dell'fft 
%      (che in Matlab corrispondono alle frequenze positive) e ottenere le restanti per simmetria.
% in realtŕ la frequenza 1 rappresenta la DC e la frequenza Nh/2+1 rappresenta la Fs e non possono essere rifasate...quindi vengono fissate a 0

%Inizializza matrici
FIR_matrix = zeros(M,No); 
%warning('OFF','MATLAB:nearlySingularMatrix')

S = No/2;    %Ritardo di causalizzazione

%Calcola R e H per ogni riga spettrale
%(Salta righe spettrali DC e Fs impostate a 0 durante l'inizializzazione)
k = 0:1:No/2;
for m = 1:M
    FIR_matrix(m,2:No/2) = (conj(IR_matrix(m,2:No/2)) .*(1/2) .* exp(-j.*2.*pi.*k(2:No/2).*S./No) ) ./ ( conj(IR_matrix(m,2:No/2)).*IR_matrix(m,2:No/2) + B(2:No/2) );  
    FIR_matrix(m,No/2+2:No) = conj( fliplr( FIR_matrix(m,2:No/2)) );
end

%Elimina matrici non piů utilizzate per risparmiare memoria
clear IR_matrix


%********************* Calcolo di fir(n) *********************

%Costruisci la finestra
twin = triang(No/2);
win = [twin(1:No/4);ones(No/2,1);twin(No/4+1:No/2)];

fir_matrix = ifft(FIR_matrix,No,2).*(ones(M,1)*win');     
clear FIR_matrix

end




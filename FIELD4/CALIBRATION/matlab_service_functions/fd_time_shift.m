%*************************************************************************************************
% fd_time_shift.m                                                          Lorenzo Chiesi 2009   *
%                                                                                                *
% Frequency Domain Time Shift                                                                                               *
% Questa funzione implementa la traslazione temporale di un segnale campionato nel dominio       *
% del tempo attraverso la moltiplicazione per un termine di  fase lineare eseguita nel dominio   *
% della frequenza.                                                                               *
% Questo metodo consente di effettuare traslazioni temporali di numeri NON interi di campioni.   * 
% Il segnale traslato mantiene la stessa lunghezza del segnale originario, la parte del segnale  *
% che fuoriesce dal vettore viene eliminata e dall'altra parte vengono introdotti degli zeri.    *       
%                                                                                                *
% out_vector = time_shift(in_vector,shift)                                                       *
%                                                                                                *
% in_vector  : Segnale campionato in ingresso                                                    *
% shift      : Numero di campioni di cui traslare il segnale                                     *
%                                                                                                *
% out_vector : Segnale traslato                                                                  *
%                                                                                                *
%*************************************************************************************************

function out_vector = fd_time_shift(in_vector,shift)

    N = length (in_vector);
    
    %Converti in un vettore colonna nel dominio della frequenza
    F_vector = fft(in_vector(:));

    %La riga spettrale 1 rappresenta la componente continua e non deve essere sfasata.
    %Applica la traslazione nel dominio della frequenza alle righe spettrali da 2 a N/2
    F_vector(2:N/2)   = F_vector(2:N/2) .* exp(-j * 2*pi * (1:1:N/2-1)' * shift/N );
    %La riga spettrale N/2+1 rappresenta la componente a frequenza Fs e poichč non puň essere sfasata viene posta a 0.    
    %F_vector(N/2+1) = 0;   
    %Le righe spettrali da N/2+2 a N hanno simmetria hermitiana con quelle da 2 a N/2
    F_vector(N/2+2:N) = conj(flipud(F_vector(2:N/2)));

    %Riconverti nel dominio del tempo
    out_vector = ifft(F_vector);
 
    %Sopprimi effetto di rotazione circolare
    z = ceil(abs(shift));  
    if (shift > 0)
        out_vector(1:z) = zeros(1,z);  
    else
        out_vector(N-z+1:N) = zeros(1,z);  
    end
    
end

%*************************************************************************************************
% spherical_directivity2.m                                                    Lorenzo Chiesi 2009 *
%                                                                                                *
% spherical_directivity(fs,sir_matrix,sir_direction,fir_matrix,theta_bar,phi_bar);                             *
%                                                                                                *
%         fs : Frequenza di campionamento                                                        *
% sir_matrix : Matrice delle IR dell'array misurate su una griglia sferica                       *
%              (dir x mic x sir_sample)   
% sir_dir : matrice delle direzioni delle risposte contenute in sir (dir x [theta_rad phi_rad]) *
% fir_matrix : Matrice dei filtri FIR dei microfoni virtuali (in x out x fir_sample)             *
% theta_bar, : Parametri opzionali che forzano la direzione in cui viene calcolata la funzione   * 
%   phi_bar    trasferimento del microfono virtuale (se omessi la direzione viene determinata    * 
%              automaticamente                                                                   *
%                                                                                                *
% Questa funzione visualizza la direttivitŕ 3D e le caratteristiche di un set di microfoni       *       
% virtuali, definiti mediante da una matrice di filtri FIR fir_matrix (in x out x fir_sample).   *          
% La funzione richiede anche l'immissione di una matrice di risposte all'impulso dell'array      *
% microfonico misurate su un piano (dir x mic x sir_sample).                                     * 
% operazioni manuali.                                                                            *
%                                                                                                *
%*************************************************************************************************


function spherical_directivity2(fs,sir_matrix,sir_dir,fir_matrix,theta_bar,phi_bar)
   

    %Calcola dimensioni della matrice delle IR dell'array microfonico
    dir        = size(sir_matrix,1);
    mic        = size(sir_matrix,2);
    sir_sample = size(sir_matrix,3);

    
    %Se viene specificata una matrice di filtri fir ne traccia la direttivitŕ,
    %altrimenti traccia la direttivitŕ delle capsule dell'array microfonico
    if exist('fir_matrix','var')
        
        
        %Calcola dimensioni della matrice dei filtri FIR
        in         = size(fir_matrix,1);
        out        = size(fir_matrix,2);
        fir_sample = size(fir_matrix,3);
        if (mic ~= in)
            error('Il numero di microfoni dell''array non corrisponde con il numero degli ingressi della matrice di filtri!')
        end

        %****************************************************************************
        % Calcola la convoluzione tra le risposte all'impulso dell'array microfonico 
        % e i filtri FIR  per ottenere la direttivitŕ dei microfoni virtuali
        %****************************************************************************
        
        %Calcola il numero di campioni necessari per avere una convoluzione NON circolare
        N = sir_sample + fir_sample ;
        
        %Trasforma i segnali in frequenza
        SIR_matrix = fft(sir_matrix,N,3);
        FIR_matrix = fft(fir_matrix,N,3);

        %Calcola la convoluzione come prodotto nel dominio della frequenza 
        %(Bastano le prime N/2+1 righe spettrali poichč i segnali sono reali)
        R_matrix = zeros(dir,out,N/2+1,'single');
        for k = 1 : N/2+1
            %Il coefficiente 2 rende la trasformata monolatera
            R_matrix(:,:,k) = 2 * SIR_matrix(:,:,k) * FIR_matrix(:,:,k);
        end
        
        %Il risultato č la matrice R_matrix (dir x out x N/2+1) che esprime 
        %per ogni microfono virtuale il guadagno in ogni direzione e ad ogni frequenza 
        
    else
    
        %Se č stata immessa solo la matrice delle risposte all'impulso dell'array microfonico
        %la matrice R_matrix (dir x mic x N/2+1) esprime per ogni microfono dell'array
        %il guadagno in ogni direzione e ad ogni frequenza
         
        out = mic;
        N = sir_sample; 
        R_matrix = fft(sir_matrix,N,3);         %Trasforma i segnali in frequenza
        R_matrix = 2*R_matrix(:,:,1:N/2+1);     %Rendi monolatera la trasformata
        
    end
  
    
    
    %****************************************************************************
    %                      Visualizzazione delle direttivitŕ 
    %****************************************************************************

    %Calcola base delle frequenze
    df = fs/N;
    f = linspace(0,fs/2,N/2+1);       

    %Calcola l'indice di alcune righe spettrali utilizzate spesso
    k_1k = round(1000/df)+1;
     
    %Crea un grafico indipendente per ogni microfono
    for out_i = 1:out

        
        %Crea finestra del grafico
        figure('Name',sprintf('Spherical Directivity (Out %d)',out_i),'NumberTitle','off')
 
        if ~exist('theta_bar','var')  
            %**********************************************************************************
            %    Calcola direzione del microfono virtuale attraverso il calcolo del baricentro
            %**********************************************************************************
            
            %rho = direttivitŕ a 1kHz
            rho1k = abs(R_matrix(:,out_i,k_1k));  
   
            %Converti in coordinate cartesiane
            [x,y,z] = sph2cart(sir_dir(:,1),sir_dir(:,2),rho1k);  

            %Calcola il valore medio secondo ogni versore (baricentro)
            xm = mean(x); 
            ym = mean(y);
            zm = mean(z);

            %Converti il baricentro in coordinate sferiche
            [theta_bar,phi_bar,r] = cart2sph(xm,ym,zm);          
        end
        
        %Calcola l'indice della direzione di caratterizzazione piů simile alla direzione della direttivitŕ massima
        [v,pos_bar] = min( (sir_dir(:,1)-theta_bar).^2 + (sir_dir(:,2)-phi_bar).^2 );


        
        %********************************************************************
        %  Risposta in frequenza e fase nella direzione del baricentro
        %********************************************************************

        %Plotta andamento del modulo della risposta in frequenza
        subplot(4,2,1)
        %semilogx(f,20*log10(abs(squeeze(R_matrix(pos_bar,out_i,:)))))
        semilogx(f,20*log10((abs(squeeze(R_matrix(pos_bar,out_i,:))))))
        ylabel('Modulo [dB]'); 
        axis([0 fs/2 -20 +10])
        set(gca,'XTick',[10,100,1000,10000,100000])
        set(gca,'XTickLAbel',{'10','100','1k','10k','100k',})
        set(gca,'XGrid','on');
        set(gca,'YGrid','on');
        title(sprintf('Modulo e fase @ 1kHz (A: %d°, E: %d°)',round(360*theta_bar/(2*pi)),round(360*phi_bar/(2*pi))))
        
        %Elimina la componente lineare dalla fase
        p = unwrap(angle(squeeze(R_matrix(pos_bar,out_i,:))));  %Estrai la risposta in fase al variare della frequenza
        i1=round(N/fs*20)+1;            %Calcola con il metodo dei minimi quadrati la componente lineare della fase (ritardo temporale)
        i2=round(N/fs*20000)+1;         %la regressione č ristretta alle frequenze da 20Hz a 20kHz che hanno un comportamento coerente
        d = f(i1:i2)'\p(i1:i2);         
        pn = p - d*f';                  %Calcola la fase a meno del ritardo
        p_1k = pn(round(1000/df)+1);    %Calcola la fase a 1kHz
        pn = pn - p_1k;                 %Poni a zero la fase f=1kHz
        pn = atan2(sin(pn),cos(pn));    %Wrap tra -pi e +pi

        %Plotta andamento della fase
        subplot(4,2,3)
        semilogx(f,pn)
        xlabel('Frequenza [Hz]');
        ylabel('Fase [rad]'); 
        axis([0 fs/2 -pi +pi])
        set(gca,'XTick',[10,100,1000,10000,100000])
        set(gca,'XTickLAbel',{'10','100','1k','10k','100k',})
        set(gca,'YTick',[-pi,-pi/2,0,pi/2,pi])
        set(gca,'YTickLAbel',{'-pi','-pi/2','0','+pi/2','+pi',})
        set(gca,'XGrid','on');
        set(gca,'YGrid','on');


        
        %********************************************************************
        %  Diagrammi di direttivitŕ sferici
        %********************************************************************

        %TODO: Ottimizzare la creazione delle matrici evitando l'uso di un ciclo
        
        %Crea le matrici nel formato necessario per la visualizzazione
        S250 = [sir_dir , zeros(dir,1) ];
        S500 = [sir_dir , zeros(dir,1) ];
        S1000 = [sir_dir , zeros(dir,1) ];
        S2000 = [sir_dir , zeros(dir,1) ];
        S4000 = [sir_dir , zeros(dir,1) ];
        S8000 = [sir_dir , zeros(dir,1) ];
        S16000 = [sir_dir , zeros(dir,1) ];

        %Scansione di tutte le direzioni
        for dir_i = 1:1:dir


                %I grafici polar3D visualizzano la fase dei numeri complessi mappandola nel colore della superfice
                %La direttivitŕ viene calcolata come il valore medio dei moduli nelle bande in 1/3 d'ottava ed č quindi un valore reale
                %La fase viene quindi sommata prendendola dalla riga spettrale di centro banda del 1/3 d'ottava

                phase =  angle( R_matrix(dir_i,out_i,round(250/df)+1)) - d*250 - p_1k;
                S250(dir_i,3)   = mean( abs( R_matrix(dir_i,out_i,round(  177/df)+1:round(  354/df)+1) ) , 3) .* exp(i*phase);

                phase =  angle( R_matrix(dir_i,out_i,round(500/df)+1)) - d*500 - p_1k;
                S500(dir_i,3)   = mean( abs( R_matrix(dir_i,out_i,round(  354/df)+1:round(  707/df)+1) ) , 3) .* exp(i*phase);

                phase =  angle( R_matrix(dir_i,out_i,round(1000/df)+1)) - d*1000 - p_1k ;
                S1000(dir_i,3)  = mean( abs( R_matrix(dir_i,out_i,round(  707/df)+1:round( 1414/df)+1) ) , 3) .* exp(i*phase);   

                phase =  angle( R_matrix(dir_i,out_i,round(2000/df)+1)) - d*2000 - p_1k;
                S2000(dir_i,3)  = mean( abs( R_matrix(dir_i,out_i,round( 1414/df)+1:round( 2828/df)+1) ) , 3) .* exp(i*phase);

                phase =  angle( R_matrix(dir_i,out_i,round(4000/df)+1)) - d*4000 - p_1k;
                S4000(dir_i,3)  = mean( abs( R_matrix(dir_i,out_i,round( 2828/df)+1:round( 5657/df)+1) ) , 3) .* exp(i*phase);

                phase =  angle( R_matrix(dir_i,out_i,round(8000/df)+1)) - d*8000 - p_1k;
                S8000(dir_i,3)  = mean( abs( R_matrix(dir_i,out_i,round( 5657/df)+1:round(11314/df)+1) ) , 3) .* exp(i*phase);

                phase =  angle( R_matrix(dir_i,out_i,round(16000/df)+1)) - d*16000 - p_1k;
                S16000(dir_i,3) = mean( abs( R_matrix(dir_i,out_i,round(11314/df)+1:N/2+1            ) ) , 3) .* exp(i*phase);
        end


        %Plotta le direttivitŕ
        AxisLimit = 1;
        
        
        subplot(4,4,9)
        magnitude = max(max(squeeze(abs(S250(:,3)))));
        S250(:,3) = S250(:,3)/magnitude;
        polar3d2(S250,AxisLimit)
        title(sprintf('250 Hz (%0.1fdB)',20*log10(magnitude)))

        subplot(4,4,10)
        magnitude = max(max(squeeze(abs(S500(:,3)))));
        S500(:,3) = S500(:,3)/magnitude;
        polar3d2(S500,AxisLimit)
        title(sprintf('500 Hz (%0.1fdB)',20*log10(magnitude)))

        subplot('position',[0.5 0.4 0.45 0.53])
        magnitude = max(max(squeeze(abs(S1000(:,3)))));
        S1000(:,3) = S1000(:,3)/magnitude;
        polar3d2(S1000,AxisLimit)
        title(sprintf('1 kHz (%0.1fdB)',20*log10(magnitude)))

        subplot(4,4,13)
        magnitude = max(max(squeeze(abs(S2000(:,3)))));
        S2000(:,3) = S2000(:,3)/magnitude;
        polar3d2(S2000,AxisLimit)
        title(sprintf('2 kHz (%0.1fdB)',20*log10(magnitude)))

        subplot(4,4,14)
        magnitude = max(max(squeeze(abs(S4000(:,3)))));
        S4000(:,3) = S4000(:,3)/magnitude;
        polar3d2(S4000,AxisLimit)
        title(sprintf('4 kHz (%0.1fdB)',20*log10(magnitude)))

        subplot(4,4,15)
        magnitude = max(max(squeeze(abs(S8000(:,3)))));
        S8000(:,3) = S8000(:,3)/magnitude;
        polar3d2(S8000,AxisLimit)
        title(sprintf('8 kHz (%0.1fdB)',20*log10(magnitude)))

        subplot(4,4,16)
        magnitude = max(max(squeeze(abs(S16000(:,3)))));
        S16000(:,3) = S16000(:,3)/magnitude;
        polar3d2(S16000,AxisLimit)
        title(sprintf('16 kHz (%0.1fdB)',20*log10(magnitude)))

    end

end


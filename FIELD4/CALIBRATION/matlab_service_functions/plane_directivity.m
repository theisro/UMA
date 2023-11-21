%*************************************************************************************************
% plane_directivity.m                                                        Lorenzo Chiesi 2009 *
%                                                                                                *
% plane_directivity(fs,hir_matrix,fir_matrix,theta_max);                                         *
%                                                                                                *
%         fs : Frequenza di campionamento                                                        *
% hir_matrix : Matrice delle IR dell'array misurate su un piano (dir x mic x hir_sample)         * 
% fir_matrix : Matrice dei filtri FIR dei microfoni virtuali (in x out x fir_sample)             *
%  theta_max : Parametro opzionale che forza la direzione in cui viene calcolata la funzione     * 
%              trasferimento del microfono virtuale (se omesso la direzione viene determinata    * 
%              automaticamente)                                                                  *
%                                                                                                *
% Questa funzione visualizza la direttivitŕ planare e le caratteristiche di un set di microfoni  *       
% virtuali, definiti mediante da una matrice di filtri FIR fir_matrix (in x out x fir_sample).   *          
% La funzione richiede anche l'immissione di una matrice di risposte all'impulso dell'array      *
% microfonico misurate su un piano (dir x mic x hir_sample).                                     * 
% operazioni manuali.                                                                            *
%                                                                                                *
%*************************************************************************************************


function plane_directivity(fs,hir_matrix,fir_matrix,theta_max)
  
    %Calcola dimensioni della matrice delle IR dell'array microfonico
    dir        = size(hir_matrix,1);
    mic        = size(hir_matrix,2);
    hir_sample = size(hir_matrix,3);

    
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
        N = hir_sample + fir_sample ;
        
        %Trasforma i segnali in frequenza
        HIR_matrix = fft(hir_matrix,N,3);
        FIR_matrix = fft(fir_matrix,N,3);

        %Calcola la convoluzione come prodotto nel dominio della frequenza 
        %(Bastano le prime N/2+1 righe spettrali poichč i segnali sono reali)
        R_matrix = zeros(dir,out,N/2+1,'single');
        for k = 1 : N/2+1
            %Il coefficiente 2 rende la trasformata monolatera
            R_matrix(:,:,k) = 2 * HIR_matrix(:,:,k) * FIR_matrix(:,:,k);
        end
        
        %Il risultato č la matrice R_matrix (dir x out x N/2+1) che esprime 
        %per ogni microfono virtuale il guadagno in ogni direzione e ad ogni frequenza 
        
    else
    
        %Se č stata immessa solo la matrice delle risposte all'impulso dell'array microfonico
        %la matrice R_matrix (dir x mic x N/2+1) esprime per ogni microfono dell'array
        %il guadagno in ogni direzione e ad ogni frequenza
        
        out = mic;
        N = hir_sample; 
        R_matrix = fft(hir_matrix,N,3);         %Trasforma i segnali in frequenza
        R_matrix = 2*R_matrix(:,:,1:N/2+1);     %Rendi monolatera la trasformata
        
    end


    
    %****************************************************************************
    %                      Visualizzazione delle direttivitŕ 
    %****************************************************************************

    
    df = fs/N;
    f = linspace(0,fs/2,N/2+1);       %Calcola base delle frequenze
    theta = (0:2*pi/dir:2*pi)';       %Calcola base degli angoli

    
    %Crea un grafico indipendente per ogni microfono
    for out_i = 1:out

        
        %Crea finestra del grafico
        figure('Name',sprintf('Planar Directivity (Out %d)',out_i),'NumberTitle','off')

        
        %************************************************
        %             Diagramma waterfall
        %************************************************
        
        %Intervallo di valori rappresentato
        dB_span = 30;
        
        %Genera mappa di colori (matrice a due dimensioni)
        map = 20*log10(squeeze(abs(R_matrix(:,out_i,:))));                %Crea una matrice quadrata (Dir x N/2+1) con i valori di direttivitŕ in dB
        map_max = max(max(map));                                    %Cerca il massimo della direttivitŕ in dB
        map = max( map , map_max-dB_span );                         %Satura il valore minimo della direttivitŕ        
        map = [map(dir/2+1:dir,:); map(1:dir/2+1,:)];               %Riorganizza l'asse degli angoli da -pi a pi 
        x_axes=[linspace(-180,0,dir/2),linspace(0,180,dir/2+1)];    %Asse degli angoli
        y_axes=linspace(0,fs/2,N/2+1);                              %Asse delle frequenze
        
        %Plotta mappa
        subplot(2,2,3)
        mesh(x_axes,y_axes,double(map'),'FaceColor','interp','FaceLighting','phong');
        axis([-180 +180 0 fs/2])
        axis square
        colorbar 
        set(gca,'Box','on');
        set(gca,'YScale','log')
        set(gca,'View',[0,90])
        set(gca,'XTick',[-180,-90,0,90,180])
        set(gca,'YTick',[16,31,63,125,250,500,1000,2000,4000,8000,16000])
        set(gca,'YTickLAbel',{'16','31','63','125','250','500','1k','2k','4k','8k','16k'})
        xlabel('Azimut [°]');
        ylabel('Frequency [Hz]');    


        

        %************************************************
        %                Diagramma polare
        %************************************************

%         %Aggrega le direttivitŕ per bande d'ottava
%         rho31 =  squeeze( mean( abs( R_matrix( : , out_i , floor(   22/df)+1:ceil(   45/df)+1 ) ) , 3) ); 
%         rho63 =  squeeze( mean( abs( R_matrix( : , out_i , floor(   45/df)+1:ceil(   89/df)+1 ) ) , 3) );
%         rho125 = squeeze( mean( abs( R_matrix( : , out_i , floor(   89/df)+1:ceil(  177/df)+1 ) ) , 3) );
%         rho250 = squeeze( mean( abs( R_matrix( : , out_i , floor(  177/df)+1:ceil(  354/df)+1 ) ) , 3) );
%         rho500 = squeeze( mean( abs( R_matrix( : , out_i , floor(  354/df)+1:ceil(  707/df)+1 ) ) , 3) );
%         rho1k  = squeeze( mean( abs( R_matrix( : , out_i , floor(  707/df)+1:ceil( 1414/df)+1 ) ) , 3) );
%         rho2k  = squeeze( mean( abs( R_matrix( : , out_i , floor( 1414/df)+1:ceil( 2828/df)+1 ) ) , 3) );
%         rho4k  = squeeze( mean( abs( R_matrix( : , out_i , floor( 2828/df)+1:ceil( 5657/df)+1 ) ) , 3) );
%         rho8k  = squeeze( mean( abs( R_matrix( : , out_i , floor( 5657/df)+1:ceil(11314/df)+1 ) ) , 3) );
%         rho16k = squeeze( mean( abs( R_matrix( : , out_i , floor(11314/df)+1:N/2+1            ) ) , 3) );
%         
        %Aggrega le direttivitŕ per bande d'ottava, Aggregazione Energetica
        rho31 =  squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor(   22/df)+1:ceil(   45/df)+1 ) ).^2 , 3) ) ); 
        rho63 =  squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor(   45/df)+1:ceil(   89/df)+1 ) ).^2 , 3) ) );
        rho125 = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor(   89/df)+1:ceil(  177/df)+1 ) ).^2 , 3) ) );
        rho250 = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor(  177/df)+1:ceil(  354/df)+1 ) ).^2 , 3) ) );
        rho500 = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor(  354/df)+1:ceil(  707/df)+1 ) ).^2 , 3) ) );
        rho1k  = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor(  707/df)+1:ceil( 1414/df)+1 ) ).^2 , 3) ) );
        rho2k  = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor( 1414/df)+1:ceil( 2828/df)+1 ) ).^2 , 3) ) );
        rho4k  = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor( 2828/df)+1:ceil( 5657/df)+1 ) ).^2 , 3) ) );
        rho8k  = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor( 5657/df)+1:ceil(11314/df)+1 ) ).^2 , 3) ) );
        rho16k = squeeze( sqrt( mean( abs( R_matrix( : , out_i , floor(11314/df)+1:N/2+1            ) ).^2 , 3) ) );
   
        %Calcola il massimo valore raggiunto dai raggi del diagramma polare
        rho_max = max([rho31; rho63; rho125; rho250; rho500; rho1k; rho2k; rho4k; rho8k; rho16k]);

        %Plotta diagramma polare normalizzando i valori
        subplot(1,2,2)
        h = polar(theta,ones(dir+1,1));               %Il primo grafico fissa la scala poi viene reso invisibile
        set(h,{'LineStyle'}, {'none'})
        hold on
        h = polar(theta,[rho31;rho31(1)]./rho_max  ,'--g');
        set(h,'LineWidth',2)
        h = polar(theta,[rho63;rho63(1)]./rho_max  ,'--b');
        set(h,'LineWidth',2)
        h = polar(theta,[rho125;rho125(1)]./rho_max,'--m'); 
        set(h,'LineWidth',2)
        h = polar(theta,[rho250;rho250(1)]./rho_max,'--y');
        set(h,'LineWidth',2)
        h = polar(theta,[rho500;rho500(1)]./rho_max,'--r');
        set(h,'LineWidth',2)
        h = polar(theta,[rho1k;rho1k(1)]./rho_max  ,'g');
        set(h,'LineWidth',2)
        h = polar(theta,[rho2k;rho2k(1)]./rho_max  ,'b');
        set(h,'LineWidth',2)
        h = polar(theta,[rho4k;rho4k(1)]./rho_max  ,'m');
        set(h,'LineWidth',2)
        h = polar(theta,[rho8k;rho8k(1)]./rho_max  ,'y');
        set(h,'LineWidth',2)
        h =polar(theta,[rho16k;rho16k(1)]./rho_max,'r');
        set(h,'LineWidth',2)       
        
        title({'- - - 31Hz green, 63Hz blue, 125Hz magenta, 250Hz yellow, 500HZ red';'------ 1kHz green, 2kHz blue, 4kHz magenta, 8kHz yellow, 16kHz red';sprintf('Gadagno a fondo scala: %0.1f dB',20*log10(rho_max))});

        
        
        
        %********************************************************************
        %  Risposta in frequenza e fase in direzione di massima direttivitŕ
        %********************************************************************
        if ~exist('theta_max','var')  
            %Calcola direzione di massima sensibilitŕ a 1 kHz
            [max_1k,pos_max] = max(rho1k);
            theta_max = (pos_max-1)*(2*pi/dir); 
        else
            pos_max = round(theta_max / (2*pi/dir)) + 1;
            max_1k = rho1k(pos_max);
        end
        %Calcola angolo a -3dB (TOFIX)
        [dummy,pos_3dB1] = min(abs(rho1k-(max_1k/sqrt(2))));
        theta_3dB1 = (pos_3dB1-1)*(2*pi/dir); 
        rho1k(pos_3dB1)=max_1k;
        [dummy,pos_3dB2] = min(abs(rho1k-(max_1k/sqrt(2))));
        theta_3dB2 = (pos_3dB2-1)*(2*pi/dir); 
        theta_3dB = sph_dist(theta_3dB1,0,theta_3dB2,0);

        %Plotta andamento del modulo della risposta in frequenza
        subplot(4,2,1)

        semilogx(f,20*log10(abs(squeeze(R_matrix(pos_max,out_i,:)))))

        ylabel('Modulo [dB]');  
        axis([0 fs/2 -20 +10])
        set(gca,'XTick',[10,100,1000,10000,100000])
        set(gca,'XTickLAbel',{'10','100','1k','10k','100k',})
        set(gca,'XGrid','on');
        set(gca,'YGrid','on');
        title(sprintf('Risposta in direzione di massima direttivitŕ (%d°) (-3dB: %d°)',round(360*theta_max/(2*pi)),round(360*theta_3dB/(2*pi))))
      
        %Elimina la componente lineare dalla fase
        p = unwrap(angle(squeeze(R_matrix(pos_max,out_i,:))));  %Estrai la risposta in fase al variare della frequenza
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


    end

end


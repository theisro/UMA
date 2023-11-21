%%********************************************************************************
%*               STORE FIR MATRIX IN A TEXT FILE FOR BRUTEFIR                    *
%*          by Lorenzo Chiesi 2009 - lorenzo.chiesi(at)gmail.com                 *
%*********************************************************************************
%
% store_brutefir(filtername,fir_matrix,gain)
%
% Input:  filtername : Name of the filter generated
%         fir_matrix : Matrix of FIR filter (In x Out x N)
%         dB_gain    : dB gain. If omitted normalization is performed

function store_brutefir(filtername,fir_matrix,dB_gain)

%Calcola parametri
in = size(fir_matrix,1);    %Numero di ingressi della matrice di filtri
out = size(fir_matrix,2);   %Numero di uscite della matrice di filtri
N = size(fir_matrix,3);     %Lunghezza dei filtri

if (in*out > 256) 
    error('BruteFir supporta fino a 256 filtri!')
end

fir_max = squeeze(max(max(max(abs(fir_matrix)))));     %Calcola il valore massimo

if not( exist('dB_gain','var') )
    
    gain = 1/fir_max;                                  %Calcola il coefficiente del guadagno  
    display(sprintf('Store BruteFir filter. (Normalization gain: %0.3f dB)', db(gain)));
else
  
    gain = 10^(dB_gain/20);
    display(sprintf('Store BruteFir filter.. (Imposed gain: %0.3f dB)',fir_file, dB_gain));
    if (abs(fir_max*gain) > 1)
        display('!!WARNING!! FILE CLIPPED!!');
    end
end



%Write coefficient file and configuration file

mkdir(filtername);

config_fid = fopen(sprintf('%s/%s_config',filtername,filtername), 'wt');

fprintf(config_fid,'#logic: "cli" { port: 3000; };\n\n');

fprintf(config_fid,'input ');
for in_i=1:in
    fprintf(config_fid,sprintf('"mic%d"',in_i));
    if (in_i < in)
        fprintf(config_fid,',');
    end
end
fprintf(config_fid,' {\n');

fprintf(config_fid,sprintf('    device: "jack" { clientname: "%s"',filtername));

%fprintf(config_fid,sprintf('    device: "jack" { clientname: "%s"; ports: ',filtername));
%for in_i=1:in
%    fprintf(config_fid,sprintf('"system:capture_%d"',in_i));
%    if (in_i < in)
%        fprintf(config_fid,',');
%    end
%end
fprintf(config_fid,'; };\n');
fprintf(config_fid,'    sample: "AUTO";\n');
fprintf(config_fid,sprintf('    channels: %d;\n',in));
fprintf(config_fid,'};\n\n');        
        
fprintf(config_fid,'output ');
for out_i=1:out
    fprintf(config_fid,sprintf('"out%d"',out_i));
    if (out_i < out)
        fprintf(config_fid,',');
    end
end
fprintf(config_fid,' {\n');

fprintf(config_fid,sprintf('    device: "jack" { clientname: "%s"',filtername));

% fprintf(config_fid,sprintf('    device: "jack" { clientname: "%s"; ports: ',filtername));
% for out_i=1:out
%     fprintf(config_fid,sprintf('"system:playback_%d"',out_i));
%     if (out_i < out)
%         fprintf(config_fid,',');
%     end
% end
fprintf(config_fid,'; };\n');
fprintf(config_fid,'    sample: "AUTO";\n');
fprintf(config_fid,sprintf('    channels: %d;\n',out));
fprintf(config_fid,'};\n\n');        



for out_i = 1:out
   for in_i = 1:in  
       
        coeff_fid = fopen(sprintf('%s/%s_%d%d',filtername,filtername,out_i,in_i), 'wt');
        for N_i = 1:N    
            fprintf(coeff_fid,'%f\n', fir_matrix(in_i,out_i,N_i));
        end
        fclose(coeff_fid);
        
        fprintf(config_fid,'coeff "spk%dmic%d" { filename: "%s_%d%d"; };\n',out_i,in_i,filtername,out_i,in_i);
        fprintf(config_fid,'filter "spk%dmic%d" {inputs: %d/6.0; outputs: %d/6.0; process: %d; coeff: "spk%dmic%d"; };\n\n',out_i,in_i,in_i-1,out_i-1,mod(out_i,2),out_i,in_i);
    end
end

fclose(config_fid);

end


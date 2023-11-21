function store_matrix(matrix_file, matrix_data)

    %Create file
    fid = fopen(matrix_file, 'w');
    
    %Write header
    fwrite(fid, length(size(matrix_data)), 'int32');
    fwrite(fid, size(matrix_data), 'int32');
    
    cells = numel(matrix_data);
    matrix_data = reshape(matrix_data,[1,cells]);
    
    %Split real and imaginary part
    matrix_data_splitted(1,:) = real(matrix_data);
    matrix_data_splitted(2,:) = imag(matrix_data);

    %Write matrix data
    fwrite(fid, matrix_data_splitted, 'single');
    
    %Close file
    fclose(fid);

end
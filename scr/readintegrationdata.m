function dataMatrix = readintegrationdata(k_value)
dataMatrix=[];
filename="integrationdata.dat";
fid = fopen(filename, 'rt'); % Apre il file in modalità lettura testo
if fid == -1
    error('Impossibile aprire il file: %s', filename);
end

line = fgetl(fid); % Legge la prima riga
found_k = false;   % Flag per indicare se la sezione 'k' è stata trovata

while ischar(line)
    % Cerca la riga 'integration degree=k'
    target_string = sprintf('integration degree=%d', k_value);
    if contains(line, target_string)
        found_k = true;
        line = fgetl(fid); % Legge la riga successiva dopo aver trovato la sezione
        
        % Importa i dati fino alla prossima sezione 'integration degree' o EOF
        current_data = [];
        while ischar(line) && ~contains(line, 'integration degree=')
            % Tenta di convertire la riga in numeri
            nums = str2num(line); %#ok<ST2NM>
            if ~isempty(nums)
                current_data = [current_data; nums]; %#ok<AGROW>
            end
            line = fgetl(fid); % Legge la riga successiva
        end
        dataMatrix = current_data; % Assegna i dati trovati
        break; % Abbiamo trovato e importato la sezione desiderata, usciamo dal loop
    else
        line = fgetl(fid); % Continua a leggere le righe
    end
end

fclose(fid); % Chiude il file

if isempty(dataMatrix) && found_k
    warning('Nessun dato numerico trovato sotto la sezione "integration degree=%d".', k_value);
elseif isempty(dataMatrix) && ~found_k
    warning('Sezione "integration degree=%d" non trovata nel file.', k_value);
end

end
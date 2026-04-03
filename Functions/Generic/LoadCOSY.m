
function [DAx] = LoadCOSY(fname,n,m,type)

% n = numero variabili dei polinomi
% m = numero di polinomi
% type = 0 ==> DA
% type = 1 ==> TM

file  = char(textread(fname,'%s','delimiter','\n'));

pos = find( (file(:,1)=='I')|(file(:,1)=='A') );
pos = pos(1:m);

DAx = [];

if n<16;
    
    for i = 1:m:length(pos)

        for j = 0:(m-1)
            k = pos(i+j)+1;
            DAx_temp(j+1).C = [];
            DAx_temp(j+1).E = [];
            while file(k,1) ~= '-'

                temp = str2num(file(k,:));
                DAx_temp(j+1).C = [ DAx_temp(j+1).C; temp(2)];
                DAx_temp(j+1).E = [ DAx_temp(j+1).E; temp(4:(4+n-1))];
                k = k+1;

            end
            if isempty(DAx_temp(j+1).C)
                DAx_temp(j+1).C = 0;
                DAx_temp(j+1).E = zeros(1,n);
            end

            if type
                
                temp = file(k+2,:);
                pos1 = find(temp=='[');
                pos2 = find(temp==',');
                pos3 = find(temp==']');
                DAx_temp(j+1).Rinf = str2double(temp(pos1+1:pos2-1));
                DAx_temp(j+1).Rsup = str2double(temp(pos2+1:pos3-1));
                
            end
            
        end

        DAx = [DAx; DAx_temp];

    end

elseif n==16
    
    for i = 1:m:length(pos)

        for j = 0:(m-1)

            k = pos(i+j)+1;
            DAx_temp(j+1).C = [];
            DAx_temp(j+1).E = [];
            while file(k,1) ~= '-'

                temp = str2num(file(k,:));
                DAx_temp(j+1).C = [ DAx_temp(j+1).C; temp(2)];
                DAx_temp(j+1).E = [ DAx_temp(j+1).E; temp(4:(4+n-1))];
                k = k+2;

            end
            if isempty(DAx_temp(j+1).C)
                DAx_temp(j+1).C = 0;
                DAx_temp(j+1).E = zeros(1,n);
            end

            if type
                
                temp = file(k+2,:);
                pos1 = find(temp=='[');
                pos2 = find(temp==',');
                pos3 = find(temp==']');
                DAx_temp(j+1).Rinf = str2double(temp(pos1+1:pos2-1));
                DAx_temp(j+1).Rsup = str2double(temp(pos2+1:pos3-1));
                
            end
            
        end

        DAx = [DAx; DAx_temp];

    end
    
else
    
    for i = 1:m:length(pos)

        for j = 0:(m-1)

            k = pos(i+j)+1;
            DAx_temp(j+1).C = [];
            DAx_temp(j+1).E = [];
            while file(k,1) ~= '-'

                temp1 = str2num(file(k,:));
                temp2 = str2num(file(k+1,:));
                DAx_temp(j+1).C = [ DAx_temp(j+1).C; temp1(2)];
                DAx_temp(j+1).E = [ DAx_temp(j+1).E; temp1(4:end-1) temp2(1:end-1)];
                
                k = k+2;

            end
            if isempty(DAx_temp(j+1).C)
                DAx_temp(j+1).C = 0;
                DAx_temp(j+1).E = zeros(1,n);
            end
            
            if type
                
                temp = file(k+2,:);
                pos1 = find(temp=='[');
                pos2 = find(temp==',');
                pos3 = find(temp==']');
                DAx_temp(j+1).Rinf = str2double(temp(pos1+1:pos2-1));
                DAx_temp(j+1).Rsup = str2double(temp(pos2+1:pos3-1));
                
            end

        end

        DAx = [DAx; DAx_temp];

    end
    
end
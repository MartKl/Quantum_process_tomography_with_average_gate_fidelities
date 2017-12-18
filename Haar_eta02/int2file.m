function int2file(int,filename)
    % writes int as text into file named filename
    text=num2str(int32(int));
    fileID = fopen(filename,'w');
    fprintf(fileID,'%s',text);
    fclose(fileID);
end

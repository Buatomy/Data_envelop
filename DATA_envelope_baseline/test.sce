Filename = "BSxxxx1--00000.csv";
File_location = "C:\Users\nut\Desktop\DATA_envelope_baseline";
for a=0:0
    if a==0  then
    else
        X_F = ascii(strsplit(Filename));
        X_F(14) = X_F(14)+1;
        Filename = ascii(X_F);
    end
    FName = fullfile(File_location,Filename);
    DATA = csvRead(FName,',',[],'string',[],[],[],5); 
    DATA = strtod(DATA);
    /*  
    T0 = find(DATA(1:$,1)>=0.05)(1);
    T1 = find(DATA(1:$,1)>=0.30)(1);*/
    XDATA_R = DATA(:,1);
    YDATA_R = DATA(:,2);
end
    plot2d(XDATA_R,YDATA_R);

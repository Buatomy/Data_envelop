clc
clear
/*
Filename = "button1--00000.csv";
File_location = "C:\Users\TUF\OneDrive - Silicon Craft Technology\URKE_verify\Sniff\RKE_Max_civic\rke_17_08_2020";
*/
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

    /*
    FName = fullfile(File_location,Filename);
    DATA = csvRead(FName,',',[],'string',[],[],[],5); 
    DATA = strtod(DATA);
    XDATA = matrix(DATA(1:$-2,1),5,-1);
    YDATA = matrix(DATA(1:$-2,2),5,-1);
    XDATA_R = XDATA(1,:);
    YDATA_R = YDATA(1,:);
    */
    //--------- Start Envelope detector --------//
    XDATA_B = XDATA_R;
    YDATA_B = YDATA_R;
    for i=1:1
        d=diff(YDATA_B)./diff(XDATA_B);
        pmax = find(d(1:$-1)>=0 & d(2:$)<0)+1;
        // edges corrections
        if d(1)<0 then
            pmax = [1 pmax]
        end
        //
        YDATA_B = YDATA_B(pmax);
        XDATA_B = XDATA_B(pmax);
    end
    X_ENV = XDATA_B;
    Y_ENV = YDATA_B;

    //disp(u)
    //--------- END Envelope detector --------//
    
    //-------- Start coumb filter (for averaging fiter) ------//
    ORD1 = 2^2;
    B = ones(1,ORD1) ;
    A = [1];
    Y_AVG = filter(B,A,Y_ENV)./ORD1;
    /// stage 2 /////
    ORD2 = 2^2;
    B = ones(1,ORD2) ;
    A = [1];
    Y_AVG = filter(B,A,Y_AVG)./ORD2;
    
    //-------- END coumb filter (for averaging fiter) ------//
    
    //-------- Start deep average for baseline ------//
    ORD3 = 2^13;
    B = ones(1,ORD3) ;
    A = [1];
    BS = filter(B,A,Y_AVG)./ORD3;
    // stage 2
    ORD4 = 2^1;
    B = ones(1,ORD4) ;
    A = [1];
    BS = filter(B,A,BS)./ORD4;
    //-------- End deep average for baseline ------//
    
    CP=0;
    //-------- Start Comparator --------- //
    for i=1:length(Y_AVG(ORD2+1:$))
        if (Y_AVG(ORD2+i)>=BS(ORD4+i)) then
            CP(1,i) = 0.14;
        else
            CP(1,i) = 0.16;
        end
    end
    //-------- End Comparator --------- //
    
    
    
    //-------- Start find Edge --------//
    //T00 = find(X_ENV>0.1)(1)
    EdgeCP = diff(CP);
    P_EdgeCP = find(EdgeCP>0);
    N_EdgeCP = find(EdgeCP<0);
    
    if length(P_EdgeCP) > length(N_EdgeCP) then
        P_EdgeCP = P_EdgeCP(1:$-1);
    elseif length(N_EdgeCP) > length(P_EdgeCP) then
        N_EdgeCP = N_EdgeCP(1:$-1);
    else ;
    end
    
    //-------- END find Edge --------//
    
    //-------- Start change time to y-axis ---------//
    
    if P_EdgeCP(1) > N_EdgeCP(1) then
        TDATA_R1 = X_ENV(P_EdgeCP)-X_ENV(N_EdgeCP);
        TDATA_R2 = X_ENV(P_EdgeCP(1:$-1))-X_ENV(N_EdgeCP(2:$));
        XR_EdgeCP = [N_EdgeCP ; P_EdgeCP];
        XRR_EdgeCP = matrix(XR_EdgeCP,1,-1);
        X_EdgeCP = [X_ENV(N_EdgeCP) ; X_ENV(P_EdgeCP)];
        X_EdgeCP = matrix(X_EdgeCP,1,-1);
        
        TDATA = [TDATA_R1 ; [TDATA_R2 ; 0]];
        TDATA1 = matrix(TDATA,1,-1);
    else
        TDATA_R1 = X_ENV(P_EdgeCP)-X_ENV(N_EdgeCP);
        TDATA_R2 = X_ENV(P_EdgeCP(1:$-1))-X_ENV(N_EdgeCP(2:$));
        XR_EdgeCP = [P_EdgeCP ; N_EdgeCP];
        XRR_EdgeCP = matrix(XR_EdgeCP,1,-1);
        X_EdgeCP = [X_ENV(P_EdgeCP) ; X_ENV(N_EdgeCP)];
        X_EdgeCP = matrix(X_EdgeCP,1,-1);
    
        TDATA = [TDATA_R1 ; [TDATA_R2 ; 0]];
        TDATA1 = matrix(TDATA,1,-1);
    end
    
    //-------- END change time to y-axis ---------//
    
    //-------- Start interpolate change time to y-axis ---------//
    DX_EdgeCP = diff(XRR_EdgeCP);
    
    for i=0:length(DX_EdgeCP)
        if i==0 then
            TDATA2 = zeros(1,XRR_EdgeCP(1)-1);
            X_Edge2 = X_ENV(1:XRR_EdgeCP(1)-1);
        else
            TRAW = ones(1,DX_EdgeCP(i)).*TDATA1(i);
            TDATA2 = [TDATA2 , TRAW];
            XRAW = X_ENV(length(X_Edge2)+1:length(X_Edge2)+DX_EdgeCP(i));
            X_Edge2 = [X_Edge2 ; XRAW];
        end
    end
    
    //-------- END intepolate point change time to y-axis ---------//
    
    
    X_Edge2 = X_Edge2';
    
    
    
    XPLearn = find((X_EdgeCP>0.1) & (X_EdgeCP<0.17) & (TDATA1>0));
    XNLearn = find((X_EdgeCP>0.1) & (X_EdgeCP<0.17) & (TDATA1<0));
    YP_Learn = median(TDATA1(XPLearn))/2;
    YN_Learn = median(TDATA1(XNLearn))/2;
    Ynstep = YN_Learn;
    Yn5pstep = YN_Learn/2;
    Ypstep = YP_Learn;
    Yp5pstep = YP_Learn/2;
    
    BDATA_T=0;
    BDATA=0;
    
    for i=1:length(TDATA1)
        if (TDATA1(i)>= Yp5pstep) && (TDATA1(i)<Ypstep+Yp5pstep) then
            BDATA_T(i) = Ypstep+0.15;
            BDATA(i) = 1;
        elseif (TDATA1(i) >= Ypstep+Yp5pstep) && (TDATA1(i) < (2*Ypstep)+Yp5pstep) then
            BDATA_T(i) = (2*Ypstep)+0.15;
            BDATA(i) = 2;
        elseif ((TDATA1(i) >= (2*Ypstep)+Yp5pstep) && (TDATA1(i) < (3*Ypstep)+Yp5pstep)) then
            BDATA_T(i) = (3*Ypstep)+0.15;
            BDATA(i) = 3;
        elseif ((TDATA1(i) >= (3*Ypstep)+Yp5pstep) && (TDATA1(i) < (4*Ypstep)+Yp5pstep)) then
            BDATA_T(i) = (4*Ypstep)+0.15;
            BDATA(i) = 4;
        elseif ((TDATA1(i) > (4*Ypstep)+Yp5pstep)) then
            BDATA_T(i) = (5*Ypstep)+0.15;
            BDATA(i) = 5;
        elseif TDATA1(i) < (4*Ynstep)+Yn5pstep then
            BDATA_T(i) = (5*Ynstep)+0.15;
            BDATA(i) = -5;
        elseif ((TDATA1(i) < (3*Ynstep)+Yn5pstep) && (TDATA1(i) > (4*Ynstep)+Yn5pstep))  then
            BDATA_T(i) = (4*Ynstep)+0.15;
            BDATA(i) = -4;
        elseif ((TDATA1(i) < (2*Ynstep)+Yn5pstep) && (TDATA1(i) > (3*Ynstep)+Yn5pstep))  then
            BDATA_T(i) = (3*Ynstep)+0.15;
            BDATA(i) = -3;
        elseif ((TDATA1(i) < Ynstep+Yn5pstep) && (TDATA1(i) > (2*Ynstep)+Yn5pstep))  then
            BDATA_T(i) = (2*Ynstep)+0.15;
            BDATA(i) = -2;
        elseif ((TDATA1(i) < Ynstep-Yn5pstep) && (TDATA1(i) > Ynstep+Yn5pstep))  then
            BDATA_T(i) = Ynstep+0.15;
            BDATA(i) = -1;
        else
            BDATA_T(i) = 0.15;
            BDATA(i) = 0;
        end
        
    end
    /*
    count = 0;
    for i=1:length(BDATA)
        
    end
    
    */
    
    
    count =0;
    for i=1:length(BDATA)
        if ((count==0) &&(BDATA(i)==2) && (BDATA(i+1)==-2) && (BDATA(i+2)==3) && (BDATA(i+3)==-4)) then
            SYC = i+4;
            count = 1;
        elseif (count==1 && ((BDATA(i)==5&&BDATA(i+1)==-5)||((BDATA(i)==-5 && BDATA(i+1)==5)))) then
            ENDMID = i-1;
            STMID = i+2;
            count = 2;
        elseif ((count==2) && (BDATA(i)==5 && BDATA(i+1)==-5)) then
            ENDlast = i-1;
            count = 3;
        elseif count == 3 then
            F_DATA = BDATA(SYC:ENDMID);
            S_DATA = BDATA(STMID:ENDlast);
            count = 4;
        else;
        end
     end
    
    for i=1:length(F_DATA)
        if i==1 then
            if F_DATA(i)>0 then
                FF_DATA = zeros(1,F_DATA(1));
            else
                FF_DATA = ones(1,-1*F_DATA(1));
            end
        else
            if F_DATA(i)>0 then
                FRAW = zeros(1,F_DATA(i));
                FF_DATA = [FF_DATA , FRAW];
            else
                FRAW = ones(1,-1*F_DATA(i));
                FF_DATA = [FF_DATA , FRAW];
            end
        end
    end
    
    for i=1:length(S_DATA)
        if i==1 then
            if S_DATA(i)>0 then
                SS_DATA = zeros(1,S_DATA(1));
            else
                SS_DATA = ones(1,-1*S_DATA(1));
            end
        else
            if S_DATA(i)>0 then
                SRAW = zeros(1,S_DATA(i));
                SS_DATA = [SS_DATA , SRAW];
            else
                SRAW = ones(1,-1*S_DATA(i));
                SS_DATA = [SS_DATA , SRAW];
            end
        end
    end
    
    clf(a)
    scf(a)
    plot2d(XDATA_R,YDATA_R);
    //plot(X_ENV,Y_ENV,"r");
    plot(X_ENV(1:$-ORD2),Y_AVG(ORD2+1:$),"r");
    plot(X_ENV(1:$-ORD4),BS(ORD4+1:$),"b");
    plot(X_ENV(1:$-ORD4-1),diff(BS(ORD4+1:$)),"g");
    plot(X_ENV(1:$-ORD2),CP,"g");
    plot(X_ENV(1:$-(ORD2+1)),EdgeCP(1:$)+0.15,"b");
    //plot(X_Edge2,TDATA2+0.15,"r");
    //plot(X_ENV(XRR_EdgeCP),BDATA_T',"m");
    //disp(u)
    //plot(X_EdgeCP,TDATA1+0.15,"r");
    plot(X_EdgeCP,TDATA1,"r")
    //plot(X_Edge2,TDATA2+0.15,"r");
    
    plot(X_Edge2,(ones(1,length(X_Edge2))*YP_Learn*4)+0.15,"b");
    plot(X_Edge2,(ones(1,length(X_Edge2))*YP_Learn*3)+0.15,"b");
    plot(X_Edge2,(ones(1,length(X_Edge2))*YP_Learn*2)+0.15,"b");
    plot(X_Edge2,(ones(1,length(X_Edge2))*YP_Learn)+0.15,"b");
    plot(X_Edge2,(ones(1,length(X_Edge2))*0.15),"black");
    plot(X_Edge2,(ones(1,length(X_Edge2))*YN_Learn)+0.15,"b");
    plot(X_Edge2,(ones(1,length(X_Edge2))*YN_Learn*2)+0.15,"b");
    plot(X_Edge2,(ones(1,length(X_Edge2))*YN_Learn*3)+0.15,"b");
    plot(X_Edge2,(ones(1,length(X_Edge2))*YN_Learn*4)+0.15,"b");
    
    
     LEN_D1(a+1) = length(FF_DATA);
     LEN_D2(a+1) = length(SS_DATA);
       
    //D1(a+1,:) = FF_DATA;
    //D2(a+1,:) = SS_DATA;
    
    
end








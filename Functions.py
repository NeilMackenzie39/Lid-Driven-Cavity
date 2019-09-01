#Functions.py
#Description: Perform discretisations for 3-step methodology
#Author: Neil Mackenzie
import numpy as np

def UCDdotC(U,C,Identifier,i,n,m):
    #Calculate Central Differenced velocity at each type of node
    if Identifier==1:                   #Bottom left corner
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=U[i,1]
        U_CD_W=U[i,0]
        U_CD_E=(U[i,0]+U[i+1,0])/2
    if Identifier ==2:                  #Bottom boundary
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=U[i,1]
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0]+U[i+1,0])/2
    if Identifier ==3:                  #Bottom right corner
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=U[i,1]
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0])     
    if Identifier ==4:                  #Left Boundary
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=U[i,0]
        U_CD_E=(U[i,0]+U[i+1,0])/2 
    if Identifier==5:                   #Internal Node
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0]+U[i+1,0])/2
    if Identifier==6:                   #Right Boundary
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0]) 
    if Identifier==7:                   #Top left corner
        U_CD_N=(U[i,1])
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0])
        U_CD_E=(U[i,0]+U[i+1,0])/2 
    if Identifier==8:                   #Top boundary
        U_CD_N=(U[i,1])
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0]+U[i+1,0])/2 
    if Identifier==9:                   #Top right corner
        U_CD_N=(U[i,1])
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0])   
        
    #Perform dot product with face normal
    U_CD_NdotC= U_CD_N*(C[i,0])
    U_CD_EdotC= U_CD_E*(C[i,1])
    U_CD_SdotC= U_CD_S*(C[i,2])
    U_CD_WdotC= U_CD_W*(C[i,3])
    
    return U_CD_NdotC,U_CD_EdotC,U_CD_SdotC,U_CD_WdotC

def NablaDotPhi(U,C,Identifier,i,n,m):
    #This function is the same as UCDotC, but excludes boundary face 
    #contributions for the convective and DeltaW* terms
    if Identifier==1:
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=0
        U_CD_W=0
        U_CD_E=(U[i,0]+U[i+1,0])/2
    if Identifier ==2:
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=0
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0]+U[i+1,0])/2
    if Identifier ==3:
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=0
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=0    
    if Identifier ==4:
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=0
        U_CD_E=(U[i,0]+U[i+1,0])/2 
    if Identifier==5:
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0]+U[i+1,0])/2
    if Identifier==6:
        U_CD_N=(U[i,1]+U[i+n,1])/2
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=0 
    if Identifier==7:
        U_CD_N=0
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=0
        U_CD_E=(U[i,0]+U[i+1,0])/2 
    if Identifier==8:
        U_CD_N=0
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=(U[i,0]+U[i+1,0])/2 
    if Identifier==9:
        U_CD_N=0
        U_CD_S=(U[i,1]+U[i-n,1])/2
        U_CD_W=(U[i,0]+U[i-1,0])/2
        U_CD_E=0
        
    U_CD_NdotC= U_CD_N*C[i,0]
    U_CD_EdotC= U_CD_E*C[i,1] 
    U_CD_SdotC= U_CD_S*C[i,2] 
    U_CD_WdotC= U_CD_W*C[i,3] 
    
    return U_CD_NdotC,U_CD_EdotC,U_CD_SdotC,U_CD_WdotC

def MatrixA(NodeEdgeLength,DeltaX,DeltaY,ControlVol,NoNodes,n,m,Identifier): 
    #This function generates Matrix A for pressure calculationusing the 
    #vertical and horizontal coefficients determined at each node
    A=np.zeros((NoNodes,NoNodes))
    for i in range (NoNodes):
        CoeffVert=NodeEdgeLength[i,1]/(ControlVol[i]*DeltaX) #East and West
        CoeffHori=NodeEdgeLength[i,0]/(ControlVol[i]*DeltaY) #North and South
       
        if Identifier[i]==1:
            #Store as 1 to use as reference pressure
            A[i,i]=1
        if Identifier[i]==2:
            A[i,i-1]=CoeffVert
            #Double vertical component due to absence of south face
            A[i,i]=-2*CoeffVert-CoeffHori
            A[i,i+1]=CoeffVert
            A[i,i+n]=CoeffHori
        if Identifier[i]==3:
            A[i,i-1]=CoeffVert
            A[i,i]=-CoeffVert-CoeffHori
            A[i,i+n]=CoeffHori
        if Identifier[i]==4:
            #Double horizontal component due to absence of west face
            A[i,i]=-CoeffVert-2*CoeffHori
            A[i,i+1]=CoeffVert
            A[i,i+n]=CoeffHori            
            A[i,i-n]=CoeffHori 
        if Identifier[i]==5:
            A[i,i-1]=CoeffVert
            A[i,i]=-2*CoeffVert-2*CoeffHori
            A[i,i+1]=CoeffVert
            A[i,i+n]=CoeffHori            
            A[i,i-n]=CoeffHori
        if Identifier[i]==6:
            A[i,i-1]=CoeffVert
            #Double horizontal component due to absence of east face
            A[i,i]=-CoeffVert-2*CoeffHori
            A[i,i+n]=CoeffHori            
            A[i,i-n]=CoeffHori
        if Identifier[i]==7:
            A[i,i]=-CoeffVert-CoeffHori
            A[i,i+1]=CoeffVert            
            A[i,i-n]=CoeffHori
        if Identifier[i]==8:
            #Double vertical component due to absence of south face
            A[i,i-1]=CoeffVert
            A[i,i]=-2*CoeffVert-CoeffHori
            A[i,i+1]=CoeffVert            
            A[i,i-n]=CoeffHori 
        if Identifier[i]==9:
            A[i,i-1]=CoeffVert
            A[i,i]=-CoeffVert-CoeffHori           
            A[i,i-n]=CoeffHori   
        
    return A     

def PDotC(P,C,Identifier,i,n,m):
    #Same as above functions, but pressure vector only has 1 component
    if Identifier==1:
        P_CD_N=(P[i]+P[i+n])/2
        P_CD_S=P[i]
        P_CD_W=P[i]
        P_CD_E=(P[i]+P[i+1])/2
    if Identifier ==2:
        P_CD_N=(P[i]+P[i+n])/2
        P_CD_S=P[i]
        P_CD_W=(P[i]+P[i-1])/2
        P_CD_E=(P[i]+P[i+1])/2
    if Identifier ==3:
        P_CD_N=(P[i]+P[i+n])/2
        P_CD_S=P[i]
        P_CD_W=(P[i]+P[i-1])/2
        P_CD_E=(P[i])    
    if Identifier ==4:
        P_CD_N=(P[i]+P[i+n])/2
        P_CD_S=(P[i]+P[i-n])/2
        P_CD_W=P[i]
        P_CD_E=(P[i]+P[i+1])/2 
    if Identifier==5:
        P_CD_N=(P[i]+P[i+n])/2
        P_CD_S=(P[i]+P[i-n])/2
        P_CD_W=(P[i]+P[i-1])/2
        P_CD_E=(P[i]+P[i+1])/2
    if Identifier==6:
        P_CD_N=(P[i]+P[i+n])/2
        P_CD_S=(P[i]+P[i-n])/2
        P_CD_W=(P[i]+P[i-1])/2
        P_CD_E=(P[i]) 
    if Identifier==7:
        P_CD_N=(P[i])
        P_CD_S=(P[i]+P[i-n])/2
        P_CD_W=(P[i])
        P_CD_E=(P[i]+P[i+1])/2 
    if Identifier==8:
        P_CD_N=(P[i])
        P_CD_S=(P[i]+P[i-n])/2
        P_CD_W=(P[i]+P[i-1])/2
        P_CD_E=(P[i]+P[i+1])/2 
    if Identifier==9:
        P_CD_N=(P[i])
        P_CD_S=(P[i]+P[i-n])/2
        P_CD_W=(P[i]+P[i-1])/2
        P_CD_E=(P[i])    
    P_CD_NdotC= P_CD_N*(C[i,0])
    P_CD_EdotC= P_CD_E*(C[i,1])
    P_CD_SdotC= P_CD_S*(C[i,2])
    P_CD_WdotC= P_CD_W*(C[i,3])
    
    return P_CD_NdotC,P_CD_EdotC,P_CD_SdotC,P_CD_WdotC
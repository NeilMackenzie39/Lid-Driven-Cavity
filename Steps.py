#Steps.py
#Description: Perform 3 steps for 2D flow solver
#Author: Neil Mackenzie
import MeshGenerator as mg
import numpy as np
import Functions as f


def Step1(U,W,NoNodes,n,m,ControlVol,Mu,DeltaX,DeltaY,Identifier,C,DeltaT):
    #Initialise matrices for step 1
    Convective = np.zeros((NoNodes,2))
    Diffusive = np.zeros((NoNodes,2))
    DeltaWStar=np.zeros((NoNodes,2))
    #Loop through each node to determine convective and diffusive term before
    #using these to determine DeltaWStar for each node
    for i in range (NoNodes):
        #Convective
        #Determine central differenced Wstar excluding boundary face contributions
        U_CD_NdotC,U_CD_EdotC,U_CD_SdotC,U_CD_WdotC=f.NablaDotPhi(U,C,Identifier[i],i,n,m)
                 
        #Find Upwind node for W    
        Upwind_N=Upwind(U_CD_NdotC,'N')
        Upwind_E=Upwind(U_CD_EdotC,'E')
        Upwind_S=Upwind(U_CD_SdotC,'S')
        Upwind_W=Upwind(U_CD_WdotC,'W')
        
        #Convective in x using W_x
        Convective[i,0]=1/ControlVol[i]*((W[i+Upwind_E,0])*(U_CD_EdotC)+(W[i+Upwind_W,0])*(U_CD_WdotC)+W[i+Upwind_S*n,0]*(U_CD_SdotC)+(W[i+Upwind_N*n,0])*(U_CD_NdotC))
        #Convective in y using W_y
        Convective[i,1]=1/ControlVol[i]*((W[i+Upwind_E,1])*(U_CD_EdotC)+(W[i+Upwind_W,1])*(U_CD_WdotC)+W[i+Upwind_S*n,1]*(U_CD_SdotC)+(W[i+Upwind_N*n,1])*(U_CD_NdotC))        
        #==============================================================================================================
        #Diffusive
        #Corners are not included here because they will always be 0
        if Identifier[i]==5:           
            #Internal Nodes include all terms
            Diffusive[i,0]=Mu*((U[i-1,0]-2*U[i,0]+U[i+1,0])/(DeltaX**2)+(U[i-n,0]-2*U[i,0]+U[i+n,0])/(DeltaY**2))
            Diffusive[i,1]=Mu*((U[i-1,1]-2*U[i,1]+U[i+1,1])/(DeltaX**2)+(U[i-n,1]-2*U[i,1]+U[i+n,1])/(DeltaY**2))
            
        if Identifier[i]==4 or Identifier[i]==6: 
            #East/West boundaries exclude west and east terms since these cancel
            Diffusive[i,0]=Mu/(DeltaY**2)*(U[i-n,0]-2*U[i,0]+U[i+n,0])
            Diffusive[i,1]=2*Mu/(DeltaX**2)*(U[i-n,1]-2*U[i,1]+U[i+n,1])
        
        if Identifier[i]==2 or Identifier[i]==8:
            #North/South boundaries exclude north and west terms since these cancel
            Diffusive[i,0]=Mu/(DeltaY**2)*(U[i-1,0]-2*U[i,0]+U[i+1,0])
            Diffusive[i,1]=2*Mu/(DeltaX**2)*(U[i+1,1]-2*U[i,1]+U[i-1,1])
    
        #Compute DeltaWStar
        DeltaWStar[i,0]=DeltaT*(Diffusive[i,0]-Convective[i,0])
        DeltaWStar[i,1]=DeltaT*(Diffusive[i,1]-Convective[i,1])        
        
    return DeltaWStar
    
def Step2(U,WStar,NoNodes,n,m,Identifier,ControlVol,C,NodeEdgeLength,DeltaX,DeltaY,rho,DeltaT,A_inv):
    #Initialise matrices for Step 2
    GradU=np.zeros((NoNodes,1))
    GradWStar=np.zeros((NoNodes,1)) 
    BMatrix=np.zeros((NoNodes,1))
    for i in range (NoNodes):
        #GradU Term
        #Determine central differenced velocities 
        U_CD_NdotC,U_CD_EdotC,U_CD_SdotC,U_CD_WdotC=f.UCDdotC(U,C,Identifier[i],i,n,m)
        GradU[i]=1/ControlVol[i]*((U_CD_EdotC)+(U_CD_WdotC)+(U_CD_NdotC)+(U_CD_SdotC))
        
        #Grad DeltaWStar Term
        #Determine central differenced Wstar excluding boundary face contributions
        WStar_CD_NdotC,WStar_CD_EdotC,WStar_CD_SdotC,WStar_CD_WdotC=f.NablaDotPhi(WStar,C,Identifier[i],i,n,m)
        GradWStar[i]=1/ControlVol[i]*((WStar_CD_EdotC)+(WStar_CD_WdotC)+(WStar_CD_NdotC)+(WStar_CD_SdotC))
        
        #Calculate B Matrix using GradU and GradDeltaWStar
        BMatrix[i]=rho/DeltaT*(GradU[i]+GradWStar[i])
    #Reset top line of B matrix for reference pressure    
    BMatrix[0]=0
    #Calculate pressure at each node using A_inverse dotted with B
    Pressure=A_inv.dot(BMatrix)
    
    return Pressure

def Step3(Pressure,Identifier,DeltaT,U,rho,n,m,C,NoNodes,ControlVol,DeltaWStar):
    #Initialise Matrices for Step 3
    U_New=np.zeros((NoNodes,2))
    GradP=np.zeros((NoNodes,2))
    #Calculate GradP in x and y-directions for each node
    for i in range(NoNodes):
        #Obtain central differenced pressures for each face
        Press_N,Press_E,Press_S,Press_W=f.PDotC(Pressure,C,Identifier[i],i,n,m)
        GradP[i,0]=1/ControlVol[i]*((Press_E)+(Press_W))
        GradP[i,1]=1/ControlVol[i]*((Press_N)+(Press_S))
    #Determine next velocity implicitly
    U_New = U - (DeltaT/rho)*GradP+(1/rho)*DeltaWStar
    
    return U_New

def Upwind(U_CDdotC,Face):
    #Determine Upwind node by checking if U dot C is positive or negative
    #where C is the area face normal
    if U_CDdotC <0:
        #Take node above or to the right if north or east face
        if Face=="N" or Face=="E":
            return 1
        #Take node below or to the left if south or west face
        if Face=="S" or Face=="W":
            return -1        
        #Use internal node for all faces if U dot C is greater than 0
    if U_CDdotC >=0:
        return 0
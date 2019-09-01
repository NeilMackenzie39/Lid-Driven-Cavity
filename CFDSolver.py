#CFDSolver.py
#Description: Solve Flow for 2D Mesh
#Author: Neil Mackenzie
import numpy as np
import MeshGenerator as mg
import Steps as st
import Graphs as gr
import Functions as f

def main(n,m):
    #Input Variables
    Mu=0.01
    rho=1
    CFL=0.5
    #Run Solver using input variables
    CFDSolver(n,m,rho,Mu,CFL)

def InitialConditions(U,rho,Identifier,NoNodes):
    #Run for loop to force velocities along top boundary (excl corners) to be 1
    #Boundary and corner velocities everywhere else are 0
    for i in range (NoNodes):
        if Identifier[i]==8:
            U[i,0]=1
            U[i,1]=0
        if Identifier[i]==1 or Identifier[i]==2 or Identifier[i]==3 or Identifier[i]==4 or Identifier[i]==6 or Identifier[i]==7 or Identifier[i]==9:
            U[i,0]=0
            U[i,1]=0
            
    #W = density X velocity at every node
    W=rho*U
    return U,W

def TimeStep(NoNodes,NodeEdgeLength,CFL,rho,Mu,U):
    #Start with large initial DeltaT
    DeltaT=10000
    #Check DeltaT at each node.
    for i in range (NoNodes):
        DeltaXSquared = (NodeEdgeLength[i,0]**2)
        DeltaYSquared = (NodeEdgeLength[i,1]**2)
        DeltaXEff=(1/(1/DeltaXSquared+1/DeltaYSquared))**0.5
        U_diff=((U[i,0]**2+U[i,1]**2)**0.5)+Mu/(rho*0.45*DeltaXEff)
        T=CFL*DeltaXEff/U_diff
    #If DeltaT at current node is smaller than current DeltaT, store new DeltaT
        if T<DeltaT:
            DeltaT=T
    return DeltaT

def Residual(U_old,U_new,Pressure_New,Pressure_Old,NoNodes):
    #Calculate residual of each velocity component as well as pressure
    #using redidual formula provided on vlua
    Resi_XVelo=((sum(((U_new[:,0])-(U_old[:,0]))**2))/NoNodes)**0.5
    Resi_YVelo=((sum(((U_new[:,1])-(U_old[:,1]))**2))/NoNodes)**0.5
    Resi_Press=((sum(((Pressure_New[:,0])-(Pressure_Old[:,0]))**2))/NoNodes)**0.5
    
    #If residual is small enough, return 'False' to terminate for loop 
    #Otherwise, return True and add 1 to no. of iterations
    if Resi_XVelo <=1*10**-5 and Resi_YVelo <=1*10**-5 and Resi_Press <=1*10**-5:
        return False,0,Resi_XVelo,Resi_YVelo,Resi_Press
    else:
        return True,1,Resi_XVelo,Resi_YVelo,Resi_Press
    
def CFDSolver(n,m,rho,Mu,CFL):
    #Run MeshGenerator to obtain Identifiers, Control Volumes, Face Normals etc
    #for each node
    Identifier,ControlVol,C,NoNodes,DeltaX,DeltaY,NodeEdgeLength=mg.NodeTable(n,m)
    
    #Initialise U and Solve initial conditions for lid-driven cavity
    U_zero=np.zeros((NoNodes,2))
    U,W=InitialConditions(U_zero,rho,Identifier,NoNodes)
    
    #Calculate Minimum Time Step
    DeltaT = TimeStep(NoNodes,NodeEdgeLength,CFL,rho,Mu,U)
    
    #Pre-Compute A Matrix for Pressures (and its inverse) to save time in while-loop
    AMatrix=f.MatrixA(NodeEdgeLength,DeltaX,DeltaY,ControlVol,NoNodes,n,m,Identifier)
    A_inv=np.linalg.inv(AMatrix)    
    
    #Start iteration and residual lists and begin with 1st iteration
    Iteration = 1
    IterationList=[]
    ResidualList=[[],[],[]]
    Run = True
    #Enter while loop to continue calculation while Run is 'True'
    #Run will be set to 'False' when the residual calculated at the end of the 
    #while loop is small enough
    while Run:
        IterationList.append(Iteration)
        
        #Run Step 1 to calculate DeltaWStar 
        DeltaWStar = st.Step1(U,W,NoNodes,n,m,ControlVol,Mu,DeltaX,DeltaY,Identifier,C,DeltaT)
        #Store current U as U_Old for residual comparison
        U_Old = U
        
        #Run Step 2
        #Set pressure = 0 if first iteration or store previous iteration pressure
        #for residual comparison
        if Iteration == 1:
            Pressure_Old=np.zeros((NoNodes,1))
        else:
            Pressure_Old = Pressure_New
        #Calculate pressure implicitly using inverse of Matrix A calculated above
        Pressure_New=st.Step2(U,DeltaWStar,NoNodes,n,m,Identifier,ControlVol,C,NodeEdgeLength,DeltaX,DeltaY,rho,DeltaT,A_inv)
        
        #Run Step 3
        #Calculate velocity at next time step
        U_New=st.Step3(Pressure_New,Identifier,DeltaT,U,rho,n,m,C,NoNodes,ControlVol,DeltaWStar)
        
        #Check Residual
        #Re-Apply Boundary Conditions
        U,W=InitialConditions(U_New,rho,Identifier,NoNodes)
        #Check residual of velocity components and pressure
        #This will return 'Run' as 'True' or 'False' depending on size of residuals
        Run,AddIteration,ResX,ResY,ResP = Residual(U_Old,U,Pressure_New,Pressure_Old,NoNodes)
        #Add 1 to iteration if it must continue, add 0 if not
        Iteration = Iteration + AddIteration
        
        #Populate residual lists for plotting
        ResidualList[0].append(ResX)
        ResidualList[1].append(ResY)
        ResidualList[2].append(ResP)
    
    #Print confirmation with number of iterations required   
    print("============================================================================================================")  
    print("Calculation Completed for " +str(m) + " by " + str (n) + " node mesh after "+ str(Iteration)+ " iterations")
        
    #Graphing Information
    U_Final=U_New
    Pressure_Final=Pressure_New[:,0]
    #Calculate magnitude using sqrt(x^2+y^2)
    U_Final_Magnitude=(U_Final[:,0]**2+U_Final[:,1]**2)**0.5
    #Colour Plots using Plotter.py
    gr.Question1Plot(U_Final_Magnitude,U_Final,Pressure_Final,n,m,IterationList,ResidualList)
    
    #Determine Velocities along x = 0.5 for Question 2
    U_X_MidPt=[]
    U_Y_MidPt=[]
    Y_Coord=[0]
    for i in range (n):
        U_X_MidPt.append(U_Final[m//2+i*n,0])
        if i!=0:
            Y_Coord.append(i/(n-1))
    for j in range (m):
        U_Y_MidPt.append(U_Final[n//2+j*m,1])
    #Plot graphs for question 2
    gr.Question2Plot(U_X_MidPt,U_Y_MidPt,Y_Coord,n,m)
 
 
#Run the main function for 21 X 21 and 51 X 51   
main(21,21)
#main(51,51)
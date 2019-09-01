#MeshGenerator.py
#Description: Generate 2D Mesh with identifiers for each node
#Author: Neil Mackenzie

import numpy as np
def NodeTable(n,m):
    
    #Define number of edges and edge length in each direction
    EdgesX=n-1  
    EdgesY=m-1
    DeltaX=1/EdgesX
    DeltaY=1/EdgesY
    #NoNodes = number nodes in x X number nodes in y
    NoNodes=n*m
    
    #Determine x and y co-ordinates along each edge
    XCoord=[0]
    YCoord=[0]
    for i in range (EdgesX):
        XCoord.append(DeltaX*(i+1))
    for j in range(EdgesY):
        YCoord.append(DeltaY*(j+1))
    
    #Generate co-ordinate list with identifier, face normal, control volume
    #and edge length for each node
    #Identifier is used to tell whether a node is a corner, edge or internal node
    #If it is a corner/edge node, Identifier shows which one
    XYCoord = np.zeros((NoNodes,2))
    Identifier=np.zeros(NoNodes)
    C=np.zeros((NoNodes,4))
    ControlVol=np.zeros(NoNodes)
    NodeEdgeLength=np.zeros((NoNodes,2))
    #Loop variable to identify location in single-row node list
    Loop=0
    #Loop through m number of nodes in x-direction to assign values to nodes
    #Then step up 1 in y-direction and repeat until y=n
    for k in range (n):
        for l in range (m):
            
            XYCoord[Loop]=XCoord[l],YCoord[k]
            if k==0 and l ==0:                              #Bottom Left Corner
                Identifier[Loop] = 1
                ControlVol[Loop] = DeltaX*DeltaY/4
                C[Loop]=DeltaX/2,DeltaY/2,-DeltaX/2,-DeltaY/2
            if k==0 and l!=0 and l!= EdgesX:                #Bottom Edge excl. corners
                Identifier[Loop] = 2 
                ControlVol[Loop] = DeltaX*DeltaY/2
                C[Loop]=DeltaX,DeltaY/2,-DeltaX,-DeltaY/2
            if k==0 and l==EdgesX:                          #Bottom Right Corner
                Identifier[Loop] = 3
                ControlVol[Loop] = DeltaX*DeltaY/4
                C[Loop]=DeltaX/2,DeltaY/2,-DeltaX/2,-DeltaY/2
            if k!=0 and l==0:                               #Left Boundary excl. corners
                Identifier[Loop] = 4
                ControlVol[Loop] = DeltaX*DeltaY/2
                C[Loop]=DeltaX/2,DeltaY,-DeltaX/2,-DeltaY
            if k!=0 and k!= EdgesY and l!= 0 and l!= EdgesX:#Internal Nodes
                Identifier[Loop] = 5
                ControlVol[Loop] = DeltaX*DeltaY
                C[Loop]=DeltaX,DeltaY,-DeltaX,-DeltaY
            if k!= 0 and k!= EdgesY and l==EdgesX:          #Right Boundary excl. corners
                Identifier[Loop] = 6 
                ControlVol[Loop] = DeltaX*DeltaY/2
                C[Loop]=DeltaX/2,DeltaY,-DeltaX/2,-DeltaY
            if k==EdgesY and l==0:                          #Top left corner
                Identifier[Loop] = 7
                ControlVol[Loop] = DeltaX*DeltaY/4
                C[Loop]=DeltaX/2,DeltaY/2,-DeltaX/2,-DeltaY/2
            if k==EdgesY and l!=0 and l!=EdgesY:            #Top boundary excl. corners
                Identifier[Loop] = 8
                ControlVol[Loop] = DeltaX*DeltaY/2
                C[Loop]=DeltaX,DeltaY/2,-DeltaX,-DeltaY/2
            if k==EdgesY and l==EdgesY:                     #Top right corner
                Identifier[Loop] = 9
                ControlVol[Loop] = DeltaX*DeltaY/4
                C[Loop]=DeltaX/2,DeltaY/2,-DeltaX/2,-DeltaY/2
                
            #Store Edge length of node. This is used to calculate time step    
            if Identifier[Loop]==1 or Identifier[Loop]==3 or Identifier[Loop]==7 or Identifier[Loop]==9:
                NodeEdgeLength[Loop]=DeltaX/2,DeltaY/2
            if Identifier[Loop]==2 or Identifier[Loop]==8:
                NodeEdgeLength[Loop]=DeltaX,DeltaY/2
            if Identifier[Loop]==4 or Identifier[Loop]==6:
                NodeEdgeLength[Loop]=DeltaX/2,DeltaY
            if Identifier[Loop]==5:
                NodeEdgeLength[Loop]=DeltaX,DeltaY
            #Add 1 to loop variable to go to next co-ordinate in list
            Loop=Loop+1
            

    return Identifier,ControlVol,C,NoNodes,DeltaX,DeltaY,NodeEdgeLength
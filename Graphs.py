#Graphs.py
#Description: Plot graphs for Question 1 and 2 of Assignment 3B
#Author: Neil Mackenzie
import Plotter as pl
import pylab as py
import matplotlib as mpl
from matplotlib.rcsetup import cycler

def Question1Plot(Velo_Magnitude,Velo,Pressure,n,m,IterationList,ResidualList):
    #Send information to Plotter.py as required by Bevan Jones' code
    Graphs=pl.Plotter("MEC4045F Assignment 3B\nPlots for a "+str(n)+" by "+str(m)+" Node Mesh",2,3,5)
    Graphs.Add1DPlot(1,"Convergence","Iteration","Residual","log",['-','-','-',],True)
    Graphs.Add2DPlot(2,"Colour Plot of Velocity Magnitude","x (m)","y (m)", m,n,True)
    Graphs.Add2DPlot(3,"Colour Plot of Pressure","x (m)","y (m)",m,n,True)
    Graphs.Add2DPlot(4,"Colour Plot of X-Velocity","x (m)","y (m)",m,n,True)
    Graphs.Add2DPlot(5,"Colour Plot of Y-Velocity","x (m)","y (m)",m,n,True)
    Legend_List = ["X-Velocity","Y-Velocity","Pressure"]
    
    Graphs.Update1DPlotData(1,IterationList,ResidualList,Legend_List)
    Graphs.Update2DPlotData(2,Velo_Magnitude,"Velocity (m/s)")
    Graphs.Update2DPlotData(3,Pressure,"Pressure (Pa)")
    Graphs.Update2DPlotData(4,Velo[:,0],"X-Velocity (m/s)")
    Graphs.Update2DPlotData(5,Velo[:,1],"Y-Velocity (m/s)")
    
def Question2Plot(U_X,U_Y,Y_Coord,n,m):
    #Set colour cycle
    mpl.rcParams['axes.prop_cycle']=cycler('color',['r','b','g'])
    #Plot X-velocity vs y-height
    py.figure(2)
    Label = str(m) +"X" + str(n) + " Node Mesh"
    py.plot(U_X,Y_Coord,label= Label)
    py.title("X-Velocity at x = 0.5")
    py.xlabel("Y-Coordinate (m)")
    py.ylabel("X-Velocity (m/s)")
    py.legend()
    
    #Plot Y-velocity vs y-height
    py.figure(3)
    Label = str(m) +"X" + str(n) + " Node Mesh"
    py.plot(U_Y,Y_Coord, label = Label)
    py.title("Y-Velocity at x = 0.5")
    py.xlabel("Y-Coordniate (m)")
    py.ylabel("Y-Velocity (m/s)")
    py.legend()
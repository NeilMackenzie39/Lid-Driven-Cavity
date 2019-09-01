# Lid-Driven-Cavity
Code developed during final year CFD course to generate a mesh and model a lid driven cavity

## Getting Started
This repository contains a jupyter notebook of the files required to execute this project

## Project Information
This project was completed for a CFD course in the final year of my mechanical engineering degree. I was one of only 3 students in the class to successfully complete the code and was awarded the class medal for the highest mark in this course. 

The inroduction below is repeated from the document submitted for this project:
The aim of this assignment was to generate code to model the case of a lid-driven cavity. This was performed using python code on equi-spaced meshes of 21 by 21 and 51 by 51 nodes. The fluid was assumed to have a density 𝜌=1 and viscosity 𝜇=0.01, and the flow domain was assumed to be a unit square. This resulted in a Reynold’s number of 𝑅𝑒=100.

The velocity was solved using the pressure projected method outlined in the course notes for MEC4045F. A fractional three-step semi-implicit method was used to solve for the velocity. First order upwinding was used for the convective term in the first step, and the pressure in the second step was solved implicitly using a coefficient matrix. This allowed for the velocity at the next time step to be solved explicitly in step three.

The discretisations for each of the terms was performed based on the method outlined in the course notes as well as that which was submitted in the author’s Assignment 3A.

## Notes on Coding
I had no experience in object-oriented programming when doing this course and therefore wrote the "functions" section of the code which identifies the location of each cell in a very inefficient manner. I knew this section of the program could be greatly improved but was running out of time for the project so I used if statements to cover each case.

## Disclaimer
The "plotter" program was written by the supervisor for the course and is the property of UCT and therefore cannot be posted here.




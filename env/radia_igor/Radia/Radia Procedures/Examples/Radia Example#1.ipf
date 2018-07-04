//=================================================
//This is the simplest example. 
//A magnetized cube is placed at position {0,0,0}. 
//It is 1 mm in size and is magnetized according to the vector  {-0.5,1,0.7} in Tesla. 
//The three components of the field at position {0.52,0.6,0.7} are computed.
//=================================================
proc RadiaExample1()
RadUtiShowHelpTopic("Example#1: Uniformly Magnetized Cube     ")

//Create necessary data structures (1D waves of 3 elements)
make/O L0 = {1,1,1}, M0 = {-0.5,1,0.7}, P0 = {0,0,0}, P1 = {0.52,0.6,0.7}, B = {0,0,0}

//Create a uniformly magnetized rectangular parallelepiped
variable m = radObjRecMag(P0,L0,M0)

//Calculate Magnetic Field at point P1
radFld(B,m,"b",P1) 
print " "
print "Magnetic field:",B[0],B[1],B[2]

KillWaves/A/Z
print " "
print "The macro generating this computation is in the file \"Radia Example#1.ipf\" accessible via menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

RadUtiShowProcExample("Radia Example#1.ipf")
RadUtiShowHelpTopic("Example#1: Uniformly Magnetized Cube     ")
end

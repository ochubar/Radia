//=================================================
//This example consists in the creation of a set of racetrack and circular coils, 
//plotting the coil geometry and the magnetic field produced. 
//This geometry corresponds to a 4T superconducting wiggler that was in operation at the ESRF. 
//It is assumed that the reader has some experience with Igor Pro and that he has successfully 
//run RadiaExample1().
//=================================================
proc RadiaExample2()
Silent 1
RadUtiShowHelpTopic("Example#2: Superconductinig Wiggler (Coils)     ")

//Current Densities in A/mm^2
variable j1=128, j2=256

//Coil Presentation Parameters
variable n1=3, n2=6, thcn=0.001
make/O c2={1,0,0}, c1={0,1,1}

//Create 5 Coils
make/O P1={0.,0.,38.}, R1={9.5,24.5}, L1={120.,0.}
variable Rt1=radObjRaceTrk(P1,R1,L1,36,n1,"man","z",j1)
radObjDrwAtr(Rt1,c1,thcn)
make/O P3={0.,0.,76.}, R3={10.,25.}, L3={90.,0.}
variable Rt3=radObjRaceTrk(P3,R3,L3,24,n1,"man","z",j1)
radObjDrwAtr(Rt3,c1,thcn)
make/O R2={24.5,55.5}, L2={120.,0.}
variable Rt2=radObjRaceTrk(P1,R2,L2,36,n1,"man","z",j2)
radObjDrwAtr(Rt2,c2,thcn)
make/O R4={25.,55.}
variable Rt4=radObjRaceTrk(P3,R4,L3,24,n1,"man","z",j2)
radObjDrwAtr(Rt4,c2,thcn)
make/O P5={0.,0.,60.}, R5={150.,166.3}, L5={0.,0.}
variable Rt5=radObjRaceTrk(P5,R5,L5,39,n2,"man","z",-j2)
radObjDrwAtr(Rt5,c2,thcn)

//Put all coils into container
make/O GrpContent={Rt1,Rt2,Rt3,Rt4,Rt5}
variable Grp=radObjCnt(GrpContent) // to remove second argument

//Define Mirror Coils
make/O P0={0,0,0}, V1={0,0,1}
radTrfZerPara(Grp,P0,V1)

//Display the geometry in a 3D viewer
//print " "
//print "ATTENTION: to continue execution of this example, please exit (close) the 3D viewer window"
radObjDrwOpenGL(Grp,"")

//Calculate and plot the Magnetic Field
make/O BzArray
make/O Pb={0,0,0}, Pe={0,300,0}
radFldLst(BzArray,Grp,"Bz",Pb,Pe,301,0)

//To calculate field, this can be also used:
//make/O/N=(3,300) Points=0; Points[1][] = q
//radFld(BzArray,Grp,"Bz",Points)

display BzArray; ModifyGraph grid(left)=1,tick(left)=2,mirror(left)=1, grid=1,tick=2,mirror=1
label bottom "Y"; label left "Bz"

//Calculate and plot the Magnetic Field Integrals
variable Npi=201, xMIn=-400, xMax=400
make/O/N=(Npi) IBzArray

SetScale/I x xMIn*0.001, xMax*0.001, "m", IBzArray
IBzArray = radFldIntCmp(Grp, "Ibz", "y", x*1000, 0, 0)

//To calculate field integral vs horizontal position, this can be used instead:
//variable x=xMIn, i=0
//variable xStep=(xMax-x)/(Npi-1)
//do
//	make/O PI1={x,-300,0}, PI2={x,300,0}, DummyIB={0}
//	IBzArray[i] = radFldInt(DummyIB,Grp,"inf","ibz",PI1,PI2)
//	i+=1; x+=xStep
//while(i<Npi)

display IBzArray; ModifyGraph grid(left)=1,tick(left)=2,mirror(left)=1, grid=1,tick=2,mirror=1
label bottom "X"; label left "Vertical Field Integral [T mm]"

TileWindows/O=1/C
KillWaves/A/Z

print " "
print "The macro generating this computation is in the file \"Radia Example#2.ipf\" accessible via menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

RadUtiShowProcExample("Radia Example#2.ipf")
RadUtiShowHelpTopic("Example#2: Superconductinig Wiggler (Coils)     ")
end


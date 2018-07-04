//=================================================
//This example creates and solves a simple Hybrid Undulator made with rectangular blocks.
//=================================================
//=================================================
//Building the Undulator
//
//The following function builds an array of rectangular permanent magnets 
//with their mirror symmetric counterparts, distributes the colors, magnetic materials 
//and sets the segmentation. The function is very generic and can be used to build 
//almost any hybrid undulator (without optimizing the design of terminations).
//=================================================
function CreateUndulator(lp, np, mp, cp, lm, nm, mm, cm, gap, gapoffset, nper)
wave lp, np, cp, lm, nm, cm
variable mp, mm, gap, gapoffset, nper

make/O zer = {0,0,0}
variable grp = radObjCnt0()

//Principal Poles and Magnets
variable y = lp[1]/4
make/O polecen = {lp[0]/4, y, -lp[2]/2 - gap/2}, poledim = {lp[0]/2, lp[1]/2, lp[2]}
variable/G pole = radObjFullMag(polecen, poledim, zer, np, grp, mp, cp)

PoleDim[1] = lp[1]
y += lp[1]/4
make/O magcen = {lm[0]/4, y, -lm[2]/2 - gap/2 - gapoffset}, magdim = {lm[0]/2, lm[1], lm[2]}
variable/G magnet = 0

variable i
for(i = 1; i <= nper; i += 1)
	make/O initm = {0, mod(i+1,2) - mod(i,2), 0}
	y += lm[1]/2
	magcen[1] = y
	magnet = radObjFullMag(magcen, magdim, initm, nm, grp, mm, cm)
	
	y += (lm[1] + lp[1])/2
	polecen[1] = y
	pole = radObjFullMag(polecen, poledim, zer, np, grp, mp, cp)
	y += lp[1]/2
endfor

make/O initm = {0, mod(nper,2) - mod(nper+1,2), 0}
y += lm[1]/4
magcen[1] = y
magdim[1] = lm[1]/2
magnet = radObjFullMag(MagCen, MagDim, InitM, nm, grp, mm, cm)

//Mirrors
make/O ex = {1,0,0}, ey = {0,1,0}, ez = {0,0,1}
radTrfZerPerp(grp, zer, ex)
radTrfZerPara(grp, zer, ez)
radTrfZerPerp(grp, zer, ey)

return grp
end

//=================================================
//Define a Magnetic Material for the Poles
//
//The following instructions define a material by means of the Magnetization  M (in Tesla) vs H 
//(in Amp/m) dependence entered point by point. The resulting M(H) curve is plotted for verification. 
//Note that the material is defined by calling the function "radMatSatIso"  with the input parameter 
//"ma" which contains the values of H and M both expressed in Tesla. The conversion between 
//Amp/m and Tesla is done by multiplying H by "4*Pi*10^(-7)" when creating "ma". 
//=================================================
function CreatePoleMaterial()
make/O H = {0.8, 1.5, 2.2, 3.6, 5, 6.8, 9.8, 18, 28, 37.5, 42, 55, 71.5, 80, 85, 88, 92, 100, 120, 150, 200, 300, 400, 600, 800, 1000, 2000, 4000, 6000, 10000, 25000, 40000}
make/O M = {0.000998995, 0.00199812, 0.00299724, 0.00499548,0.00699372, 0.00999145, 0.0149877, 0.0299774, 0.0499648, 0.0799529, 0.0999472, 0.199931, 0.49991, 0.799899, 0.999893, 1.09989, 1.19988, 1.29987, 1.41985, 1.49981, 1.59975, 1.72962, 1.7995, 1.89925, 1.96899, 1.99874,  2.09749, 2.19497, 2.24246, 2.27743, 2.28958, 2.28973}
make/O/N=(2, dimsize(H, 0)) ma
ma[0] = (4*Pi*1.e-7)*H[q]
ma[1] = M[q]
return radMatSatIso(ma, "tab")
end

//=================================================
//Convenience function to facilitate plotting of M vs H curves
//=================================================
function RadAuxMatMvsH(obj, cmp, hx, hy, hz)
variable obj, hx, hy, hz
string cmp
make/O h = {hx, hy, hz}, m = {0, 0, 0}
return radMatMvsH(m, obj, cmp, h)
end

//=================================================
//Create the Geometry, Solve the Problem and Calculate the Magnetic Field
//
//The following instructions build an undulator with NdFeB magnets and poles made with the material 
//previously defined. One could also use a Radia pre-defined steel by replacing "mp = CreatePoleMaterial()" 
//with "mp = radMatStd("AFK502", 0)" or "mp = radMatStd("Xc06", 0)" or "mp = radMatStd("Steel42", 0)". 
//After the materials are created, the function building the undulator structure is called. 
//Then the geometry is solved in order to determine magnetization vectors in all magnets and poles;
//after this, the magnetic field is computed.
//=================================================
proc RadiaExample3()
Silent 1
RadUtiShowHelpTopic("Example#3: Generic Hybrid Undulator     ")

//General Parameters
variable gap = 20, nper = 2, per = 46, gapoffset = 1

//Pole Parameters
make/O lp = {45,5,25}, np = {2,2,5}, cp = {1,0,1}
variable ll = per/2 - lp[1]

//Magnet Parameters
make/O lm = {65,ll,45}, nm = {1,3,1}, cm = {0,1,1}

//Cleaning the Memory
radUtiDelAll()

//Defining the Pole Material
variable mp = CreatePoleMaterial()
 
//Defining the Permanent Magnet Material: NdFeB with 1.2 Tesla Remanent Magnetization
variable mm = radMatStd("NdFeB", 1.2)

//Building the Structure
variable grp = CreateUndulator(lp, np, mp, cp, lm, nm, mm, cm, gap, gapoffset, nper)

//Displaying the Geometry in a 3D Viewer
radObjDrwOpenGL(grp,"")
//radObjDrwQD3D(grp,"")

//Solving the Problem for Magnetization
make/O/N=4 res
make/O zer = {0,0,0}, B = {0,0,0}

variable t0 = startmstimer
radSolve(res, grp, 0.0003, 1000, 4)

variable t1 = stopmstimer(t0)/1e6
print "Solved in :    ", round(t1 - t0), "  seconds"
print "Number of Iterations : ", res[3]
print "Average Stability of Magnetization at the last iteration :  ", res[0], "  T"
print "Maximum Absolute Magnetization at the last iteration : ", res[1], "  T"
print "Maximum H vector at the last iteration : ", res[2], "  T"
print ""
print "Central Field  Bz(0,0,0) : ", radFld(B, grp, "Bz", zer), "  T"

//Calculating the Magnetic Field
variable y1 = -(nper+2)/2*per, y2 = (nper+2)/2*per, npt = 101
variable sy = (y2 - y1)/(npt - 1)
make/O/N=(3, npt) pts; pts = 0
pts[1][] = y1 + sy*q
make/O bz
radFld(bz, grp, "bz", pts)

//Plotting the Magnetic Field
setscale/I x y1*0.001, y2*0.001, "m", bz
display bz
modifygraph grid(left)=1,tick(left)=2,mirror(left)=1, grid=1,tick=2,mirror=1
label bottom "y"; label left "Bz"
textbox/C/N=text0/A=LT "On-Axis Magnetic Field"

//Checking the Magnetic Material Characteristics
variable nh = 101, h1 = -0.01, h2 = 0.01
variable sh =  (h2 - h1)/(nh - 1)
make/O/N=(nh) mzp
mzp = RadAuxMatMvsH(pole, "mz", 0, 0, h1 + sh*p)
setscale/I x h1, h2, "T", mzp
setscale/I d -2, 2, "T", mzp
display mzp
modifygraph grid(left)=1,tick(left)=2,mirror(left)=1, grid=1,tick=2,mirror=1
label bottom "H [T]"; label left "M [T]"
textbox/C/N=text0/A=LT "Pole Material Characteristics"

nh = 101; h1 = -1; h2 = 1
sh =  (h2 - h1)/(nh - 1)
make/O/N=(nh) mzm, mym
mym = RadAuxMatMvsH(magnet, "my", 0, h1 + sh*p, 0)
mzm = RadAuxMatMvsH(magnet, "mz", 0, 0, h1 + sh*p)
setscale/I x h1, h2, "T", mym, mzm
setscale/I d -2, 2, "T", mym, mzm
display mym, mzm
modifygraph grid(left)=1,tick(left)=2,mirror(left)=1, grid=1,tick=2,mirror=1
label bottom "H [T]"; label left "M [T]"
textbox/C/N=text0/A=LT "NdFeB Characteristics"

TileWindows/O=1/C
//KillWaves/A/Z

print " "
print "The macro generating this computation is in the file \"Radia Example#3.ipf\" accessible via menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

RadUtiShowProcExample("Radia Example#3.ipf")
RadUtiShowHelpTopic("Example#3: Generic Hybrid Undulator     ")
end

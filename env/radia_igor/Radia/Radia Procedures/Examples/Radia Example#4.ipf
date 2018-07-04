//=================================================
//This example illustrates the use of a polyhedron shape by means of the "radObjMltExtPgn" function. 
//A uniformly magnetized polyhedron is created with a magnetization of 1 Tesla. The field produced by 
//this polyhedron is computed and shown to be uniform inside the volume of the polyhedron and equal to 
//2/3 Tesla as expected from an analytical integration.
//=================================================
//=================================================
//Function Building the Geometry
//
//These instructions define a uniformly magnetized convex polyhedron which converges to a sphere at large 
//values of segmentation parameters.
//=================================================
function CreateMagnetizedSphere(r, nphi, nz, m)
variable r, nphi, nz
wave m

make/O/D slices = {{0, 0}}
make/O/N=(nz + 1) numpts; numpts[0] = 1
make/O/N=(nz + 1) attitudes; attitudes[0] = -r

variable dz = 2.*r/nz
variable z = -r + dz, dphi = 2.*Pi/nphi
variable i, k, theta, costheta, phi
variable pcount = 1, pcountpr = 0

for(i = 1; i < nz; i += 1)
	theta = asin(z/r); costheta = cos(theta); phi = dphi
	pcountpr = pcount; pcount += 1
	redimension/D/N=(2, pcount) slices; slices[0][pcountpr] = r*costheta; slices[1][pcountpr] = 0
	for(k = 1; k < nphi; k += 1)
		pcountpr = pcount; pcount += 1
		redimension/D/N=(2, pcount) slices; slices[0][pcountpr] = r*cos(phi)*costheta; slices[1][pcountpr] = r*sin(phi)*costheta
		phi += dphi
	endfor
	numpts[i] = nphi
	attitudes[i] = z
	z += dz
endfor
pcountpr = pcount; pcount += 1
redimension/N=(2, pcount) slices; slices[0][pcountpr] = 0; slices[1][pcountpr] = 0
numpts[nz] = 1
attitudes[nz] = r

return radObjMltExtPgn(slices, numpts, attitudes, m)
end

//=================================================
//Create the Geometry, Calculate and Plot the Magnetic Field
//=================================================
proc RadiaExample4()
Silent 1
//RadUtiShowHelpTopic("Example#4: Uniformly Magnetized Polyhedron     ")

radUtiDelAll()
make/O m = {1.,0.,0.}, color = {0,0.5,0.8}
variable sphere = CreateMagnetizedSphere(1., 15, 15, m)
radObjDrwAtr(sphere, color, 0.001)

//Display the geometry in a 3D viewer
//print " "
//print "ATTENTION: to continue execution of this example, please exit (close) the 3D viewer window"
radObjDrwOpenGL(sphere, "")

//Calculating the Magnetic Field
make/O pt = {0,0,0}, b = {0,0,0}, bx
radFld(b, sphere, "b", pt)
print "Field in the Center [T] : ", b[0], b[1], b[2]

variable y1 = -1.1, y2 = 1.1, npt = 101
variable sy = (y2 - y1)/(npt - 1)
make/O/N=(3, npt) pts; pts = 0; pts[1][] = y1 + sy*q
radFld(bx, sphere, "bx", pts)

//Plotting the Magnetic Field
setscale/I x y1*0.001, y2*0.001, "m", bx
display bx
modifygraph grid(left)=1,tick(left)=2,mirror(left)=1, grid=1,tick=2,mirror=1
label bottom "y"; label left "Bx"
textbox/C/N=text0/A=LT "Magnetic Field Inside the Sphere"

TileWindows/O=1/C
//KillWaves/A/Z

print " "
print "The macro generating this computation is in the file \"Radia Example#4.ipf\" accessible via menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

RadUtiShowProcExample("Radia Example#4.ipf")
RadUtiShowHelpTopic("Example#4: Uniformly Magnetized Polyhedron     ")
end

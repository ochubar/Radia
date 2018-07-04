

//:Evaluate:      BeginPackage["Radia`"]


//-------------------------------------------------------------------------


:Evaluate:      radObjRecMag::usage = "radObjRecMag[{x,y,z},{wx,wy,wz},{mx,my,mz}:{0,0,0}] creates a rectangular parallelepipedic block with center point {x,y,z}, dimensions {wx,wy,wz}, and magnetization {mx,my,mz}."
:Evaluate:      radObjThckPgn::usage = "radObjThckPgn[x,lx,{{y1,z1},{y2,z2},...},a:\"x\",{mx,my,mz}:{0,0,0}] creates an extruded polygon block; x is the position of the block's center of gravity in the extrusion direction, lx is the thickness, {{y1,z1},{y2,z2},...} is a list of points describing the polygon in 2D; the extrusion direction is defined by the character a (which can be \"x\", \"y\" or \"z\"), {mx,my,mz} is the block magnetization."
:Evaluate:      radObjPolyhdr::usage = "radObjPolyhdr[{{x1,y1,z1},{x2,y2,z2},...},{{f1i1,f1i2,...},{f2i1,f2i2,...},...},{mx,my,mz}:{0,0,0},J->{jx,jy,jz}|{{jx,jy,jz},{{djxdy,djxdy,djxdz},{djydy,djydy,djydz},{djzdy,djzdy,djzdz}}},Lin->Rel|Abs] creates a polyhedron. {{x1,y1,z1},{x2,y2,z2},...} is a list of the polyhedron vertex points, {{f1n1,f1n2,...},{f2n1,f2n2,...},...} is a list of lists of indexes of vertex points defining the polyhedron faces, {mx,my,mz} is magnetization inside the polyhedron. The option J->... can be used to define constant {jx,jy,jz} or linearly-varying current density vector inside the polyhedron; the linear dependence can be defined through 3x3 matrix of coefficients {{djxdy,djxdy,djxdz},{djydy,djydy,djydz},{djzdy,djzdy,djzdz}}; depending on the value of the option Lin->... this linear dependence is treated with respect to the polyhedron center (Lin->Rel, default) or with respect to the origin of the Cartesian frame (Lin->Abs)."
:Evaluate:      radObjMltExtPgn::usage = "radObjMltExtPgn[{{{{x11,y11},{x12,y12},...},z1},{{{x21,y21},{x22,y22},...},z2},...},{mx,my,mz}:{0,0,0}] attempts to create a convex polyhedron or a set of convex polyhedrons based on horizontal slices, which should be convex planar polygons. The slice polygons are defined by the nested list {{{{x11,y11},{x12,y12},...},z1},{{{x21,y21},{x22,y22},...},z2},...}, with {{x11,y11},{x12,y12},...},... describing the polygons in 2D, and z1, z2,... giving their attitudes (vertical coordinates). {mx,my,mz} is the magnetization inside the polyhedron(s) created."
:Evaluate:      radObjMltExtPgnCur::usage = "radObjMltExtPgnCur[z:0,a:\"z\",{{{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...}},I,Frame->Loc|Lab] attempts to create a set of current-carrying convex polyhedron objects by applying a generalized extrusion to the initial planar convex polygon. The initial polygon is defined for the \"attitude\" z (0 by default) by the list of 2D points {{x1,y1},{x2,y2},...}, with the  a  character specifying orientation of this polygon normal in 3D space: if a = \"z\" (default orientation), the polygon is assumed to be parallel to XY plane of the laboratory frame (\"y\" for ZX plane, \"x\" for YZ plane). The extrusion can consist of a number of \"steps\", with each step creating one convex polyhedron defined optionally by one (or combination of) rotation(s), and/or translation(s), and/or one homothety: {Rk,Tk,Hk}, k = 1,2,..., applied to the base polygon (i.e. either the initial base polygon, or the polygon obtained by previous extrusion step). In case if k-th extrusion step contains one rotation Rk about an axis, it is defined as {{xRk,yRk,zRk},{vxRk,vyRk,vzRk},phRk}}, where {xRk,yRk,zRk} and {vxRk,vyRk,vzRk} are respectively 3D coordinates of a point and a vector difining the rotation axis, and phRk the rotation angle in radians; in case if Rk is a combination of \"atomic\" rotations about different axes, it should be defined as list: {Rk1,Rk2,...}. If k-th extrusion step includes translation Tk, it must be defined by vector {vxTk,vyTk,vzTk}; optional homothety with respect to the base polygon center of gravity should be defined either by two different coefficients {pxHk,pyHk} with respect to two orthogonal axes of the base polygon local frame, or by nested list {{pxHk,pyHk},phHk}, where phHk is rotation angle of the two homothety axes in radians. A real number I defines average current in Amperes along the extrusion path. The Frame->Loc or Frame->Lab option specifies whether the transformations at each step of the extrusion path are defined in the frame of the previous base polygon (Frame->Loc, default), or all the transformations are defined in the laboratory frame (Frame->Lab)."
:Evaluate:      radObjMltExtPgnMag::usage = "radObjMltExtPgnMag[z:0,a:\"z\",{{{x1,y1},{x2,y2},...},{{k1,q1},{k2,q2},...}:{{1,1},{1,1},...},{{R1,T1,H1},{R2,T2,H2},...}},{{mx1,my1,mz1},{mx2,my2,mz2},...}:{{0,0,0},{0,0,0},...},Frame->Loc|Lab,ki->Numb|Size,TriAngMin->...,TriAreaMax->...,TriExtOpt->\"...\"] attempts to create a set of uniformly magnetized polyhedron objects by applying a generalized extrusion to the initial planar convex polygon. The initial polygon is defined for the \"attitude\" z (0 by default) by the list of 2D points {{x1,y1},{x2,y2},...}, with the  a  character specifying orientation of this polygon normal in 3D space: if a = \"z\" (default orientation), the polygon is assumed to be parallel to XY plane of the laboratory frame (\"y\" for ZX plane, \"x\" for YZ plane). The extrusion can consist of a number of \"steps\", with each step creating one convex polyhedron defined optionally by one (or combination of) rotation(s), and/or translation(s), and/or one homothety: {Rk,Tk,Hk}, k = 1,2,..., applied to the base polygon (i.e. either the initial base polygon, or the polygon obtained by previous extrusion step). In case if k-th extrusion step contains one rotation Rk about an axis, it is defined as {{xRk,yRk,zRk},{vxRk,vyRk,vzRk},phRk}}, where {xRk,yRk,zRk} and {vxRk,vyRk,vzRk} are respectively 3D coordinates of a point and a vector difining the rotation axis, and phRk the rotation angle in radians; in case if Rk is a combination of \"atomic\" rotations about different axes, it should be defined as a list: {Rk1,Rk2,...}. If k-th extrusion step includes translation Tk, it must be defined by vector {vxTk,vyTk,vzTk}; optional homothety with respect to the base polygon center of gravity should be defined either by two different coefficients {pxHk,pyHk} with respect to two orthogonal axes of the base polygon local frame, or by nested list {{pxHk,pyHk},phHk}, where phHk is rotation angle of the two homothety axes in radians. Optional list {{mx1,my1,mz1},{mx2,my2,mz2},...} defines magnetization vectors in each of polyhedrons to be created. The Frame->Loc or Frame->Lab option specifies whether the transformations at each step of the extrusion path are defined in the frame of the previous base polygon (Frame->Loc, default), or all the transformations are defined in the laboratory frame (Frame->Lab). Optionally, the object can be subdivided by (extruded) triangulation at its creation; this occurs if {{k1,q1},{k2,q2},...} subdivision (triangulation) parameters for each segment of the base polygon border are defined; the meaning of k1, k2,... depends on value of the option ki: if ki->Numb (default), then k1, k2,... are subdivision numbers; if ki->Size, they are average sizes of sub-segments to be produced; q1, q2,... are ratios of the last-to-first sub-segment lengths; the TriAngMin option defines minimal angle of triangles to be produced (in degrees, default is 20); the TriAreaMax option defines maximal area of traingles to be produced (in mm^2, not defined by default); the ExtOpt option allows to specify additional parameters for triangulation function in a string."
:Evaluate:      radObjMltExtRtg::usage = "radObjMltExtRtg[{{{x1,y1,z1},{wx1,wy1}},{{x2,y2,z2},{wx2,wy2}},...},{mx,my,mz}:{0,0,0}] attempts to create a convex polyhedron or a set of convex polyhedrons based on horizontal slices of rectangular shape. The slice rectangles are defined by the nested list {{{x1,y1,z1},{wx1,wy1}},{{x2,y2,z2},{wx2,wy2}},...}, with {x1,y1,z1}, {x2,y2,z2},... being center points of the rectangles, and {wx1,wy1}, {wx2,wy2},... their dimensions. {mx,my,mz} is the magnetization inside the polyhedron(s) created."
:Evaluate:      radObjMltExtTri::usage = "radObjMltExtTri[x,lx,{{y1,z1},{y2,z2},...},{{k1,q1},{k2,q2},...},a:\"x\",{mx,my,mz}:{0,0,0},ki->Numb|Size,TriAngMin->...,TriAreaMax->...,TriExtOpt->\"...\"] creates triangulated extruded polygon block, i.e. a straight prism with its bases being (possibly non-convex) polygons subdivided by triangulation; x is the position of the block's center of gravity in the extrusion direction, lx is the thickness, {{y1,z1},{y2,z2},...} is a list of points describing the polygon in 2D; {{k1,q1},{k2,q2},...} are subdivision (triangulation) parameters for each segment of the base polygon border; the meaning of k1, k2,... depends on value of the option ki: if ki->Numb (default), then k1, k2,... are subdivision numbers; if ki->Size, they are average sizes of sub-segments to be produced; q1, q2,... are ratios of the last-to-first sub-segment lengths; the extrusion direction is defined by the character a (which can be \"x\", \"y\" or \"z\"); {mx,my,mz} is magnetization inside the block. The TriAngMin option defines minimal angle of triangles to be produced (in degrees, default is 20); the TriAreaMax option defines maximal area of traingles to be produced (in mm^2, not defined by default); the ExtOpt option allows to specify additional parameters for triangulation function in a string."
:Evaluate:      radObjArcPgnMag::usage = "radObjArcPgnMag[{x,y},a,{{r1,z1},{r2,z2},...},{phimin,phimax},nseg,\"sym|nosym\":\"nosym\",{mx,my,mz}:{0,0,0}] creates a finite-length arc of polygonal cross-section with the position and orientation of the rotation axis defined by pair of coordinates {x,y} and character a (which can be \"x\", \"y\" or \"z\"), the cross-section 2D polygon {{r1,z1},{r2,z2},...}, initial and final rotation angles {phimin,phimax}, number of sectors (segments) vs azimuth nseg, and magnetization vector {mx,my,mz}. Depending on the value of the \"sym|nosym\" switch, the magnetization vectors in nseg sector polyhedrons are either assumed to satisfy rotational symmetry conditions (\"sym\"), or are assumed independent (i.e. will be allowed to vary independently at further relaxation)."
:Evaluate:      radObjCylMag::usage = "radObjCylMag[{x,y,z},r,h,nseg,a:\"z\",{mx,my,mz}:{0,0,0}] creates a cylindrical magnet approximated by a right polygon with center point {x,y,z}, base radius r, height h, number of segments nseg, orientation of the rotation axis defined by character a (which can be \"x\", \"y\" or \"z\"), and magnetization vector {mx,my,mz}."
:Evaluate:      radObjRecCur::usage = "radObjRecCur[{x,y,z},{wx,wy,wz},{jx,jy,jz}] creates a current carrying rectangular parallelepipedic block with center point {x,y,z}, dimensions {wx,wy,wz}, and current density {jx,jy,jz}."
:Evaluate:      radObjArcCur::usage = "radObjArcCur[{x,y,z},{rmin,rmax},{phimin,phimax},h,nseg,j,\"man|auto\":\"man\",a:\"z\"] creates a current-carrying finite-length arc of rectangular cross-section, center point {x,y,z}, inner and outer radii {rmin,rmax}, initial and final angles {phimin,phimax}, height h, number of segments nseg, and azimuthal current density j. According to the value of the \"man|auto\" switch, the field from the arc is computed based on the number of segments nseg (\"man\"), or on the general absolute precision level specified by the function radFldCmpCrt (\"auto\"). The orientation of the rotation axis is defined by the character a (which can be either \"x\", \"y\" or \"z\")."
:Evaluate:      radObjRaceTrk::usage = "radObjRaceTrk[{x,y,z},{rmin,rmax},{lx,ly},h,nseg,j,\"man|auto\":\"man\",a:\"z\"] creates a current carrying racetrack coil consisting of four 90-degree bents connected by four straight parts of rectangular straight section, center  point {x,y,z}, inner and outer bent radii {rmin,rmax}, straight section lengths {lx,ly}, height h, number of segments in bents nseg, and azimuthal current density j. According to the value of the \"man|auto\" switch, the field from the bents is computed based on the number of segments nseg (\"man\"), or on the general absolute precision level specified by the function radFldCmpCrt (\"auto\"). The orientation of the racetrack axis is defined by the character a (which can be either \"x\", \"y\" or \"z\")."
:Evaluate:      radObjFlmCur::usage = "radObjFlmCur[{{x1,y1,z1},{x2,y2,z2},...},i] creates a filament polygonal line conductor defined by the sequence of points {{x1,y1,z1},{x2,y2,z2},...} with current i."
:Evaluate:      radObjCnt::usage = "radObjCnt[{obj1,obj2,...}] creates a container object for the objects {obj1,obj2,...}."
:Evaluate:      radObjAddToCnt::usage = "radObjAddToCnt[cnt,{obj1,obj2,...}] adds objects {obj1,obj2,...} to the container object cnt."
:Evaluate:      radObjCntStuf::usage = "radObjCntStuf[obj] gives a list of general indexes of the objects present in container if obj is a container; or returns {obj} if obj is not a container."
:Evaluate:      radObjBckg::usage = "radObjBckg[{bx,by,bz}] creates a source of uniform background magnetic field {bx,by,bz}."
:Evaluate:      radObjDivMag::usage = "radObjDivMag[obj,{{k1,q1},{k2,q2},{k3,q3}},{\"pln\",{n1x,n1y,n1z},{n2x,n2y,n2z},{n3x,n3y,n3z}}|{\"cyl\",{{ax,ay,az},{vx,vy,vz}},{px,py,pz},rat},kxkykz->Numb|Size,Frame->Loc|Lab|LabTot] subdivides the object obj. The main subdivision parameters are defined by the list {{k1,q1},{k2,q2},{k3,q3}}. The meaning of k1, k2 and k3 depends on the value of the option kxkykz: if kxkykz->Numb (default), then k1, k2 and k3 are subdivision numbers; if kxkykz->Size, they are average sizes of the sub-objects to be produced. q1, q2 and q3 in any case are ratios of the last-to-first sub-object sizes. The third variable is optional. If it is not used, the subdivision is performed in directions X, Y and Z. If it is used in the form {\"pln\",...}, the function performs subdivision by three sets of parallel planes normal to the vectors {n1x,n1y,n1z}, {n2x,n2y,n2z} and {n3x,n3y,n3z}. The distances between the parallel planes are defined by the parameters {k1,q1},{k2,q2} and {k3,q3}. If the third variable is used in the form {\"cyl\",...}, the function performs subdivision by a system of coaxial elliptic cylinders. The cylinder axis is defined by the point {ax,ay,az} and vector {vx,vy,vz}. One of two axes of the cylinder base ellipses is exactly the perpendicular from the point {px,py,pz} to the cylinder axis; rat is the ratio of the ellipse axes lengths. In the case of the subdivision by elliptic cylinders, the parameters {k1,q1},{k2,q2} and {k3,q3} correspond to radial, azimuthal, and axial directions respectively. If the option Frame is set to Frame->Loc (default), the subdivision is performed in local frame(s) of the object(s); if Frame->Lab or Frame->LabTot, the subdivision is performed in the laboratory frame. The actions of Frame->Lab and Frame->LabTot differ for containers only: Frame->Lab means that each of the objects in the container is subdivided separately; Frame->LabTot means that the objects in the container are subdivided as one object, by the same planes."
:Evaluate:      radObjCutMag::usage = "radObjCutMag[obj,{x,y,z},{nx,ny,nz},Frame->Loc|Lab|LabTot] cuts the object obj by a plane normal to the vector {nx,ny,nz} and passing through the point {x,y,z}. The Frame->Loc, Frame->Lab or Frame->LabTot option specifies whether the cuting plane is defined in the local frame of the object obj or in the laboratory frame (default). The actions of Frame->Lab and Frame->LabTot differ for containers only: Frame->Lab means that each of the objects in the container is cut separately; Frame->LabTot means that the objects in the container are cut as one object, by the same plane. The function returns a list of indexes of the objects produced by the cutting."
:Evaluate:      radObjGeoVol::usage = "radObjGeoVol[obj] computes geometrical volume of the object obj."
:Evaluate:      radObjGeoLim::usage = "radObjGeoLim[obj] computes geometrical limits of the object obj in the laboratory frame. Returns {xmin, xmax, ymin, ymax, zmin, zmax}."
:Evaluate:      radObjDpl::usage = "radObjDpl[obj,FreeSym->False|True] duplicates the object obj. The option FreeSym->False|True specifies whether the symmetries (transformations with multiplicity more than one) previously applied to the object obj should be simply copied at the duplication (FreeSym->False, default), or a container of new independent objects should be created in place of any symmetry previously applied to the object obj. In both cases the final object created by the duplication has exactly the same geometry as the initial object obj."
:Evaluate:      radObjDrwAtr::usage = "radObjDrwAtr[obj,{r,g,b},thcn] applies RGB color {r,g,b} and line thickness thcn to object obj."
:Evaluate:      radObjDrw::usage = "radObjDrw[obj] prepares a set of 3D graphical primitives representing object obj in space."
//:Evaluate:      radObjDrwQD3D::usage = "radObjDrwQD3D[obj,EdgeLines->True|False,Faces->True|False,Axes->True|False] starts an application for viewing 3D geometry of the object obj in interactive mode. The feature is implemented using QuickDraw 3D graphics library from Apple Computer. The option EdgeLines->True|False (default EdgeLines->True) highlights the edge lines of objects; the option Faces->True|False (default Faces->True) shows faces of the objects; the option Axes->True|False (default Axes->True) shows the Cartesian frame axes."
:Evaluate:      radObjDrwOpenGL::usage = "radObjDrwOpenGL[obj,EdgeLines->True|False,Faces->True|False,Axes->True|False] starts an application for viewing 3D geometry of the object obj. The viewer is based on the GLUT / OpenGL graphics library. The option EdgeLines->True|False (default EdgeLines->True) highlights the edge lines of objects; the option Faces->True|False (default Faces->True) shows faces of the objects; the option Axes->True|False (default Axes->True) shows the Cartesian frame axes."
:Evaluate:      radObjDegFre::usage = "radObjDegFre[obj] gives number of degrees of freedom for the relaxation of the object obj."
:Evaluate:      radObjM::usage = "radObjM[obj] returns coordinates of geometrical center point and magnetization vector of the object obj at that point. If obj is a container, a list of the container members' center points and their magnetizations is returned."
:Evaluate:      radObjCenFld::usage = "radObjCenFld[obj,\"A|B|H|J|M\"] returns coordinates of geometrical center point of the object obj and a field characteristic vector at that point. The type of field characteristic is defined by the second parameter (character); it can be one of the following: \"A\" for vector potential, \"B\" for magnetic field induction, \"H\" for magnetic field strength, \"J\" for current density, \"M\" for magnetization. If obj is a container, a list of the container members' center points and their field characteristics is returned."
:Evaluate:      radObjScaleCur::usage = "radObjScaleCur[obj,k] scales current (density) in the obj by multiplying it by k (if obj is a current-carying object). If obj is a container, the current (density) scaling applies to all its members."
:Evaluate:      radTrfTrsl::usage = "radTrfTrsl[{vx,vy,vz}] creates a translation with vector {vx,vy,vz}."
:Evaluate:      radTrfRot::usage = "radTrfRot[{x,y,z},{vx,vy,vz},phi] creates a rotation of angle phi around the axis defined by the point {x,y,z} and the vector {vx,vy,vz}."
:Evaluate:      radTrfPlSym::usage = "radTrfPlSym[{x,y,z},{nx,ny,nz}] creates a symmetry with respect to plane defined by point {x,y,z} and normal vector {nx,ny,nz}."
:Evaluate:      radTrfInv::usage = "radTrfInv[] creates a field inversion."
:Evaluate:      radTrfCmbL::usage = "radTrfCmbL[OrigTrf,trf] multiplies original transformation OrigTrf by another transformation trf from left."
:Evaluate:      radTrfCmbR::usage = "radTrfCmbR[OrigTrf,trf] multiplies original transformation OrigTrf by another transformation trf from right."
:Evaluate:      radTrfOrnt::usage = "radTrfOrnt[obj,trf] orients object obj by applying transformation trf to it once."
:Evaluate:      radTrfMlt::usage = "radTrfMlt[obj,trf,mlt] creates mlt-1 objects. Each object is derived from the previous one by applying the transformation trf. Following this, the object obj becomes equivalent to mlt different objects."
:Evaluate:      radMatLin::usage = "radMatLin[{ksipar,ksiper},mr] or radMatLin[{ksipar,ksiper},{mrx,mry,mrz}] creates a linear anisotropic magnetic material with susceptibilities parallel (perpendicular) to the easy magnetization axis given by ksipar (ksiper). In the first form of the function, mr is the magnitude of the remanent magnetization vector; the direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is applied (the magnetization vector is specified at the object creation). In the second form, {mrx,mry,mrz} is the remanent magnetization vector explicitly defining the direction of the easy magnetization axis in any object to which the material is later applied."
:Evaluate:      radMatSatIso::usage = "radMatSatIso[{ksi1,ms1},{ksi2,ms2},{ksi3,ms3}] creates a nonlinear isotropic magnetic material with the magnetization magnitude M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3) where H is the magnitude of the magnetic field strength vector. radMatSatIso[{{H1,M1},{H2,M2},...}] creates a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs {{H1,M1},{H2,M2},...} in Tesla."
:Evaluate:      radMatSatLam::usage = "radMatSatLam[data,p,{nx,ny,nz}:0] creates laminated nonlinear anisotropic magnetic material with stacking factor p and vector normal to the lamination planes given by {nx,ny,nz}. data is a list of pairs of numbers defining the material dependence M(H) as for isotropic case. If this list contains only 3 elements, it is interpreted as {{ksi1,ms1},{ksi2,ms2},{ksi3,ms3}} and M(H) is then computed by the formula M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3); if it contains more than 3 elements, it is interpreted as tabulated dependence M(H): {{H1,M1},{H2,M2},...} in Tesla. If the vector normal {nx,ny,nz} is not supplied, the lamination planes are assumed to be perpendicular to the magnetization vector in the object to which the material is applied (the magnetization vector is specified at the object creation)."
:Evaluate:      radMatSatAniso::usage = "radMatSatAniso[datapar,dataper] where datapar can be {{ksi1,ms1,hc1},{ksi2,ms2,hc2},{ksi3,ms3,hc3},{ksi0,hc0}} or ksi0, and dataper can be {{ksi1,ms1},{ksi2,ms2},{ksi3,ms3},ksi0} or ksi0 - creates a nonlinear anisotropic magnetic material. If the first argument is set to {{ksi1,ms1,hc1},{ksi2,ms2,hc2},{ksi3,ms3,hc3},{ksi0,hc0}}, the magnetization vector component parallel to the easy axis is computed as ms1*tanh(ksi1*(hpa-hc1)/ms1) + ms2*tanh(ksi2*(hpa-hc2)/ms2) + ms3*tanh(ksi3*(hpa-hc3)/ms3) + ksi0*(hpa-hc0), where hpa is the field strength vector component parallel to the easy axis. If the second argument is set to {{ksi1,ms1},{ksi2,ms2},{ksi3,ms3},ksi0}, the magnetization vector component perpendicular to the easy axis is computed as ms1*tanh(ksi1*hpe/ms1) + ms2*tanh(ksi2*hpe/ms2) + ms3*tanh(ksi3*hpe/ms3) + ksi0*hpe, where hpe is the field strength vector component perpendicular to the easy axis. If the first or second argument is set to ksi0, the magnetization component parallel (perpendicular) to the easy axis is computed by ksi0*hp, where hp is the corresponding component of field strength vector. At least one of the magnetization vector components should non-linearly depend on the field strength. The direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is later applied."
:Evaluate:      radMatApl::usage = "radMatApl[obj,mat] applies material mat to object obj."
:Evaluate:      radMatMvsH::usage = "radMatMvsH[obj,\"mx|my|mz\"|\"\",{hx,hy,hz}] computes magnetization from magnetic field strength vector (hx,hy,hz) for the material of the object obj; the magnetization components are specified by the second argument."
:Evaluate:      radRlxPre::usage = "radRlxPre[obj,srcobj:0] builds an interaction matrix for the object obj, treating the object srcobj as additional external field source."
:Evaluate:      radRlxMan::usage = "radRlxMan[intrc,meth,iternum,rlxpar] executes manual relaxation procedure for interaction matrix intrc using method number meth, by making iternum iterations with relaxation parameter value rlxpar."
:Evaluate:      radRlxAuto::usage = "radRlxAuto[intrc,prec,maxiter,meth:4,ZeroM->True|False] executes automatic relaxation procedure with the interaction matrix intrc using the method number meth. Relaxation stops whenever the change in magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than maxiter. The option value ZeroM->True (default) starts the relaxation by setting the magnetization values in all paricipating objects to zero; ZeroM->False starts the relaxation with the existing magnetization values in the sub-volumes."
:Evaluate:      radRlxUpdSrc::usage = "radRlxUpdSrc[intrc] updates external field data for the relaxation (to take into account e.g. modification of currents in coils, if any) without rebuilding the interaction matrix."
:Evaluate:      radSolve::usage = "radSolve[obj,prec,maxiter,meth:4] builds interaction matrix for the object obj and executes relaxation procedure using the method number meth. Relaxation stops whenever the change in magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than maxiter."
:Evaluate:      radFldCmpCrt::usage = "radFldCmpCrt[prcB,prcA,prcBint,prcFrc,prcTrjCrd,prcTrjAng] sets general absolute accuracy levels for computation of field induction (prcB), vector potential (prcA), induction integrals along straight line (prcBint), field force (prcFrc), relativistic particle trajectory coordinates (prcTrjCrd) and angles (prcTrjAng)."
:Evaluate:      radFldCmpPrc::usage = "radFldCmpPrc[PrcB->prb,PrcA->pra,PrcBInt->prbint,PrcForce->prfrc,PrcTorque->prtrq,PrcEnergy->pre,PrcCoord->prcrd,PrcAngle->prang] sets general absolute accuracy levels for computation of magnetic field induction, vector potential, induction integral along straight line, field force, torque, energy; relativistic charged particle trajectory coordinates and angles. The function works in line with the Mathematica mechanism of Options. PrcB, PrcA, PrcBInt, PrcForce, PrcTorque, PrcEnergy, PrcCoord, PrcAngle are names of the options; prb, pra, prbint, prfrc, prtrq, pre, prcrd, prang are the corresponding values (real numbers specifying the accuracy levels)."
:Evaluate:      radFld::usage = "radFld[obj,\"bx|by|bz|hx|hy|hz|ax|ay|az|mx|my|mz\"|\"\",{x,y,z}|{{x1,y1,z1},{x2,y2,z2},...}] computes magnetic field created by the object obj in point(s) {x,y,z} ({x1,y1,z1},{x2,y2,z2},...). The field component is specified by the second input variable. The function accepts a list of 3D points of arbitrary nestness: in this case it returns the corresponding list of magnetic field values."
:Evaluate:      radFldLst::usage = "radFldLst[obj,\"bx|by|bz|hx|hy|hz|ax|ay|az|mx|my|mz\"|\"\",{x1,y1,z1},{x2,y2,z2},np,\"arg|noarg\":\"noarg\",strt:0.] computes magnetic field created by object obj in np equidistant points along a line segment from {x1,y1,z1} to {x2,y2,z2}; the field component is specified by the second input variable; the \"arg|noarg\" string variable specifies whether to output a longitudinal position for each point where the field is computed, and strt gives the start-value for the longitudinal position."
:Evaluate:      radFldInt::usage = "radFldInt[obj,\"inf|fin\",\"ibx|iby|ibz\"|\"\",{x1,y1,z1},{x2,y2,z2}] computes magnetic field induction integral produced by the object obj along a straight line specified by points {x1,y1,z1} and {x2,y2,z2}; depending on the second variable value, the integral is infinite (\"inf\") or finite from {x1,y1,z1} to {x2,y2,z2} (\"fin\"); the field integral component is specified by the third input variable. The unit is Tesla x millimeter."
:Evaluate:      radFldFrc::usage = "radFldFrc[obj,shape] computes force of the field produced by the object obj into a shape defined by shape. shape can be the result of radObjRecMag[..] (parallelepiped) or radFldFrcShpRtg[..] (rectangular surface)."
:Evaluate:      radFldFrcShpRtg::usage = "radFldFrcShpRtg[{x,y,z},{wx,wy}] creates a rectangle with central point {x,y,z} and dimensions {wx,wy}; to be used for field force computation."
:Evaluate:      radFldEnr::usage = "radFldEnr[objdst,objsrc] or radFldEnr[objdst,objsrc,{kx,ky,kz}] computes potential energy (in Joule) of the object objdst in the field created by the object objsrc. The first form of the function performes the computation based on absolute accuracy value for the energy (by default 10 Joule; can be modified by the function radFldCmpPrc). The second form performs the computation based on the destination object subdivision numbers {kx,ky,kz}."
:Evaluate:      radFldEnrFrc::usage = "radFldEnrFrc[objdst,objsrc,\"fx|fy|fz|\"\"] or radFldEnrFrc[objdst,objsrc,\"fx|fy|fz|\"\",{kx,ky,kz}] computes force (in Newton) acting on the object objdst in the field produced by the object objsrc. The first form of the function performes the computation based on absolute accuracy value for the force (by default 10 Newton; can be modified by the function radFldCmpPrc). The second form performs the computation based on the destination object subdivision numbers {kx,ky,kz}."
:Evaluate:      radFldEnrTrq::usage = "radFldEnrTrq[objdst,objsrc,\"tx|ty|tz|\"\",{x,y,z}] or radFldEnrTrq[objdst,objsrc,\"tx|ty|tz|\"\",{x,y,z},{kx,ky,kz}] computes torque (in Newton*mm) with respect to point {x,y,z}, acting on the object objdst in the field produced by the object objsrc. The first form of the function performes the computation based on absolute accuracy value for the torque (by default 10 Newton*mm; can be modified by the function radFldCmpPrc). The second form performs the computation based on the destination object subdivision numbers {kx,ky,kz}."
:Evaluate:      radFldPtcTrj::usage = "radFldPtcTrj[obj,E,{x0,dxdy0,z0,dzdy0},{y0,y1},np] computes transverse coordinates and its derivatives (angles) of a relativistic charged particle trajectory in 3D magnetic field produced by the object obj, using the 4th order Runge-Kutta integration. The particle energy is E [GeV], initial transverse coordinates and derivatives are {x0,dxdy0,z0,dzdy0}; the longitudinal coordinate y is varied from y0 to y1 in np steps. All positions are in millimeters and angles in radians."
:Evaluate:      radFldFocPot::usage = "radFldFocPot[obj,{x1,y1,z1},{x2,y2,z2},np] computes the potential for trajectory of relativistic charged particle in magnetic field produced by the object obj. The integration is made from {x1,y1,z1} to {x2,y2,z2} with np equidistant points."
:Evaluate:      radFldFocKickPer::usage = "radFldFocKickPer[obj,{x1,y1,z1},{nsx,nsy,nsz},per,nper,{n1x,n1y,n1z},r1,np1,r2,np2,com:\"\",{nh:1,nps:8,d1:0,d2:0},\"T2m2|rad|microrad\":\"T2m2\",en:1,\"fix|tab\":\"fix\"] computes matrices of 2nd order kicks of trajectory of relativistic charged particle in periodic magnetic field produced by the object obj. The longitudinal integration along one period starts at point {x1,y1,z1} and is done along direction pointed by vector {nsx,nsy,nsz}; per is period length, nper is number of full periods; one direction of the transverse grid is pointed by vector {n1x,n1y,n1z}, the other transverse direction is given by vector product of {n1x,n1y,n1z} and {nsx,nsy,nsz}; r1 and r2 are ranges of the transverse grid, np1 and np2 are corresponding numbers of points; com is arbitrary string comment; nh is maximum number of magnetic field harmonics to treat (default 1), nps is number of longitudinal points (default 8), d1 and d2 are steps of transverse differentiation (by default equal to the steps of the transverse grid); the \"T2m2|rad|microrad\" string variable specifies the units for the resulting 2nd order kick values (default \"T2m2\"); en is electron elergy in GeV (optional, required only if units are \"rad\" or \"microrad\"); the \"fix|tab\":\"fix\" string variable specifies the format of the output data string, \"fix\" for fixed-width, \"tab\" for tab-delimited (i.e. element [[6]] of the output list, default \"fix\"). Returns list containing: [[1]]- matrix of kick values in the first transverse direction, [[2]]- matrix of kick values in the second transverse direction, [[3]]- matrix of longitudinally-integrated squared transverse magnetic field calculated on same transverse mesh as kicks, [[4]],[[5]]- lists of positions defining the transverse grid, [[6]]- formatted string containing the computed results (for saving into a text file)."
:Evaluate:      radFldFocKick::usage = "radFldFocKick[obj,{x1,y1,z1},{nsx,nsy,nsz},{ds1,ds2,ds3,...},nps,{n1x,n1y,n1z},r1,np1,r2,np2,com:\"\",{d1:0,d2:0}] computes matrices of 2nd order kicks of trajectory of relativistic charged particle in arbitrary magnetic field produced by the object obj. PLEASE NOTE that this is a time-consuming calculation; therefore try using the radFldFocKickPer function were applicable. The longitudinal integration starts at point {x1,y1,z1} and is done along direction pointed by vector {nsx,nsy,nsz}; the transverse matrices of the kick values are computed for the planes located at distances {ds1,ds2,ds3,...} from {x1,y1,z1}; nps is total number of points for longitudinal integration; one direction of the transverse grid is pointed by vector {n1x,n1y,n1z}, the other transverse direction is given by vector product of {n1x,n1y,n1z} and {nsx,nsy,nsz}; r1 and r2 are ranges of the transverse grid, np1 and np2 are corresponding numbers of points; com is arbitrary string comment; d1 and d2 are steps of transverse differentiation (by default equal to the steps of the transverse grid). Returns list containing: [[1]]- list of triplets of 2D matrices representing kick values in the first and second transverse directions and the longitudinally-integrated squared transverse magnetic field, including the total matrices computed for whole range of longitudinal position and partial matrices corresponding to given longitudinal positions, [[2]],[[3]]- lists of positions defining the transverse grid, [[4]]- formatted string containing all the computed results (for saving into a text file)."
:Evaluate:      radFldShimSig::usage = "radFldShimSig[obj,\"bx|by|bz|hx|hy|hz|ibx|iby|ibz\"|\"\",{dx,dy,dz},{x1,y1,z1},{x2,y2,z2},np,{vix,viy,viz}:{0,0,0}] computes virtual \"shim signature\", i.e. variation of magnetic field component defined by the second variable, introduced by displacement {dx,dy,dz} of magnetic field source object obj. The field variation is computed at np equidistant points along a line segment from {x1,y1,z1} to {x2,y2,z2}; the vector {vix,viy,viz} is taken into account if a field integral variation should be computed: in this case, it defines orientation of the integration line." 
:Evaluate:      radFldLenTol::usage = "radFldLenTol[abs,rel,zero:0] sets absolute and relative randomization magnitudes for all the length values, including coordinates and dimensions of the objects producing magnetic field, and coordinates of points where the field is computed. Optimal values of the variables can be: rel=10^(-11), abs=L*rel, zero=abs, where L is the distance scale value (in mm) for the problem to be solved. Too small randomization magnitudes can result in run-time code errors."
:Evaluate:      radFldLenRndSw::usage = "radFldLenRndSw[\"on|off\"] switches on or off the randomization of all the length values. The randomization magnitude can be set by the function radFldLenTol."
:Evaluate:      radFldUnits::usage = "radFldUnits[] shows the physical units currently in use."
:Evaluate:      radUtiDel::usage = "radUtiDel[elem] deletes element elem."
:Evaluate:      radUtiDelAll::usage = "radUtiDelAll[] deletes all previously created elements."
:Evaluate:      radUtiDmp::usage = "radUtiDmp[elem,\"asc|bin\":\"asc\"] outputs information about elem, which can be either one element or list of elements; second argument specifies whether the output should be in ASCII (\"asc\", default) or in Binary (\"bin\") format."
:Evaluate:      radUtiDmpPrs::usage = "radUtiDmpPrs[bstr] parses byte-string bstr produced by radUtiDmp[elem,\"bin\"] and attempts to instantiate elem objects(s); returns either index of one instantiated object (if elem was an index of one object) or a list of indexes of instantiated objects (if elem was a list of objects)."
:Evaluate:      radUtiIntrptTim::usage = "radUtiIntrptTim[t] sets interruption time quanta in seconds for platforms with no preemptive multitasking."
:Evaluate:      radUtiVer::usage = "radUtiVer[] returns the version number of the Radia executable file."


:Evaluate:      rAdObjPgn::usage = "rAdObjPgn[z,{{x1,y1},{x2,y2},...},{mx,my,mz}:{0,0,0}] "
:Evaluate:      rAdObjSetLocMgn::usage = "rAdObjSetLocMgn[obj,{{{ix1,iy1,iz1},{mx1,my1,mz1}},{{ix2,iy2,iz2},{mx2,my2,mz2}},...}] "
:Evaluate:      rAdObjDrwWithoutTrfMlt::usage = "rAdObjDrwWithoutTrfMlt[obj] "
:Evaluate:      rAdObjDelDrwAtr::usage = "rAdObjDelDrwAtr[obj] "
:Evaluate:      rAdObjGeoVol::usage = "rAdObjGeoVol[obj] "
:Evaluate:      rAdObjRecMagsAsExtPgns::usage = "rAdObjRecMagsAsExtPgns[\"on|off\"] "
:Evaluate:      rAdObjRecMagsAsPolyhdrs::usage = "rAdObjRecMagsAsPolyhdrs[\"on|off\"] "
:Evaluate:      rAdObjExtPgnsAsPolyhdrs::usage = "rAdObjExtPgnsAsPolyhdrs[\"on|off\"] "
:Evaluate:      rAdObjRcgnRecMags::usage = "rAdObjRcgnRecMags[\"on|off\"] "
:Evaluate:      rAdObjDivByPlns::usage = "rAdObjDivByPlns[obj,{{{n1x,n1y,n1z},k1,q1},{{n2x,n2y,n2z},k2,q2},...},\"lab|loc\":\"lab\"] "
:Evaluate:      rAdObjDplFrSym::usage = "rAdObjDplFrSym[obj] duplicates the geometry of the object obj in such a way that a container of new independent objects appears instead of any symmetry (i.e., transformation with multiplicity more than one) previously applied to the object obj."
:Evaluate:      rAdRlxOutIntrcMatr::usage = "rAdRlxOutIntrcMatr[intrc] "
:Evaluate:      rAdRlxOutIntrcVect::usage = "rAdRlxOutIntrcVect[intrc,\"ext|tot|mag\"] "
:Evaluate:      rAdFldMltpolThr::usage = "rAdFldMltpolThr[{a0,a1,a2,a3}] "
:Evaluate:      rAdFldCmpMetNxNyNz::usage = "rAdFldCmpMetNxNyNz[obj,1|0] "
:Evaluate:      rAdUtiDelAll::usage = "rAdUtiDelAll[] "
:Evaluate:      rAdUtiRlxMemAllocMet::usage = "rAdUtiRlxMemAllocMet[\"tot|parts\"] "
:Evaluate:      rAdUtiRetInp::usage = "rAdUtiRetInp[inp,ntimes] "
:Evaluate:      rAdUtiStrtProf::usage = "rAdUtiStrtProf[Flag,NumFuncs,Depth] "
:Evaluate:      rAdUtiStopProf::usage = "rAdUtiStopProf "


//:Evaluate:      Begin["Radia`Private`"]


//-------------------------------------------------------------------------


void RecMag P(( double,double,double, double,double,double, double,double,double ));

:Begin:
:Function:      RecMag
:Pattern:       radObjRecMag[{xCoordin_,yCoordin_,zCoordin_}, {wxWidth_,wyWidth_,wzWidth_}, Magnetiz_List:{0.,0.,0.}]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[wxWidth],N[wyWidth],N[wzWidth], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
:ArgumentTypes: { Real,Real,Real, Real,Real,Real, Real,Real,Real }
:ReturnType:    Manual
:End:


void ExtrudedPolygon P(( ));

:Begin:
:Function:      ExtrudedPolygon
:Pattern:       radObjThckPgn[xCoordin_, lxWidth_, ListOf2dPoints_, Magnetiz_List:{0.,0.,0.}, Orient_String:"x"]
:Arguments:     { N[xCoordin], N[lxWidth], N[ListOf2dPoints], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]], Orient }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void ExtrudedPolygon2 P(( ));

:Begin:
:Function:      ExtrudedPolygon2
:Pattern:       radObjThckPgn[xCoordin_, lxWidth_, ListOf2dPoints_, Orient_String:"x", Magnetiz_List:{0.,0.,0.}]
:Arguments:     { N[xCoordin], N[lxWidth], N[ListOf2dPoints], Orient, N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//:Begin:
//:Function:      ExtrudedPolygon2
//:Pattern:       radObjThckPgnMag[xCoordin_, lxWidth_, ListOf2dPoints_, Orient_String:"x", Magnetiz_List:{0.,0.,0.}]
//:Arguments:     { N[xCoordin], N[lxWidth], N[ListOf2dPoints], Orient, N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void PlanarPolygon P(( ));

:Begin:
:Function:      PlanarPolygon
:Pattern:       rAdObjPgn[zCoordin_, ListOf2dPoints_, Magnetiz_List:{0.,0.,0.}]
:Arguments:     { N[zCoordin], N[ListOf2dPoints], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void Polyhedron1 P(( ));

:Begin:
:Function:      Polyhedron1
:Pattern:       radObjPolyhdr[ListOfPoints_, ListOfListOfIndexes_, OptPar1_:0, OptPar2_:0, OptPar3_:0]
:Arguments:     { N[ListOfPoints], Round[ListOfListOfIndexes], N[OptPar1], N[OptPar2], N[OptPar3] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:

//:Begin:
//:Function:      Polyhedron1
//:Pattern:       radObjPolyhdr[ListOfPoints_, ListOfListOfIndexes_, Magnetiz_List:{0,0,0}, OptPar1_:0]
//:Arguments:     { N[ListOfPoints], Round[ListOfListOfIndexes], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]], N[OptPar1_] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:

//:Begin:
//:Function:      Polyhedron1
//:Pattern:       radObjPolyhdrMag[ListOfPoints_, ListOfListOfIndexes_, Magnetiz_List:{0,0,0}]
//:Arguments:     { N[ListOfPoints], Round[ListOfListOfIndexes], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void RecMagsAsExtrPolygons P(( const char* ));

:Begin:
:Function:      RecMagsAsExtrPolygons
:Pattern:       rAdObjRecMagsAsExtPgns[OnOrOff_String]
:Arguments:     { OnOrOff }
:ArgumentTypes: { String }
:ReturnType:    Manual
:End:


void RecMagsAsPolyhedrons P(( const char* ));

:Begin:
:Function:      RecMagsAsPolyhedrons
:Pattern:       rAdObjRecMagsAsPolyhdrs[OnOrOff_String]
:Arguments:     { OnOrOff }
:ArgumentTypes: { String }
:ReturnType:    Manual
:End:


void ExtPgnsAsPolyhedrons P(( const char* ));

:Begin:
:Function:      ExtPgnsAsPolyhedrons
:Pattern:       rAdObjExtPgnsAsPolyhdrs[OnOrOff_String]
:Arguments:     { OnOrOff }
:ArgumentTypes: { String }
:ReturnType:    Manual
:End:


void RecognizeRecMags P(( const char* ));

:Begin:
:Function:      RecognizeRecMags
:Pattern:       rAdObjRcgnRecMags[OnOrOff_String]
:Arguments:     { OnOrOff }
:ArgumentTypes: { String }
:ReturnType:    Manual
:End:


void MultGenExtrPolygon P(( ));

:Begin:
:Function:      MultGenExtrPolygon
:Pattern:       radObjMltExtPgn[ListOfLayerPolygons_, Magnetiz_List:{0,0,0}]
:Arguments:     { N[ListOfLayerPolygons], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//:Begin:
//:Function:      MultGenExtrPolygon
//:Pattern:       radObjMltExtPgnMag[ListOfLayerPolygons_, Magnetiz_List:{0,0,0}]
//:Arguments:     { N[ListOfLayerPolygons], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void MultGenExtrPolygonCur P(( ));
//radObjMltExtPgnCur[z:0,a:"z",{{x1,y1},{x2,y2},...},{{R1,T1,H1},{R2,T2,H2},...},I]

:Begin:
:Function:      MultGenExtrPolygonCur
:Pattern:       radObjMltExtPgnCur[zCoordin_, Orient_, ListOf2dPoints_:0, ExtrusionPath_:0, AvgCurrent_:0, OptPar1_:0]
:Arguments:     { N[zCoordin], N[Orient], N[ListOf2dPoints], N[ExtrusionPath], N[AvgCurrent], N[OptPar1] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void MultGenExtrPolygonMag P(( ));
//radObjMltExtPgnMag[z:0,a:"z",{{x1,y1},{x2,y2},...},{{k1,q1},{k2,q2},...},{{R1,T1,H1},{R2,T2,H2},...},{{mx1,my1,mz1},{mx2,my2,mz2},...}]

:Begin:
:Function:      MultGenExtrPolygonMag
:Pattern:       radObjMltExtPgnMag[zCoordin_, Orient_, ListOf2dPoints_:0, ListOfSubdivParam_:0, ExtrusionPath_:0, Magnetiz_:0, OptPar1_:0, OptPar2_:0, OptPar3_:0, OptPar4_:0, OptPar5_:0]
:Arguments:     { N[zCoordin], N[Orient], N[ListOf2dPoints], N[ListOfSubdivParam], N[ExtrusionPath], N[Magnetiz], OptPar1, OptPar2, OptPar3, OptPar4, OptPar5  }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void MultGenExtrRectangle P(( ));

:Begin:
:Function:      MultGenExtrRectangle
:Pattern:       radObjMltExtRtg[ListOfLayerRectangles_, Magnetiz_List:{0,0,0}]
:Arguments:     { N[ListOfLayerRectangles], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//:Begin:
//:Function:      MultGenExtrRectangle
//:Pattern:       radObjMltExtRtgMag[ListOfLayerRectangles_, Magnetiz_List:{0,0,0}]
//:Arguments:     { N[ListOfLayerRectangles], N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void MultGenExtrTriangle P(( ));

:Begin:
:Function:      MultGenExtrTriangle
:Pattern:       radObjMltExtTri[xCoordin_, lxWidth_, ListOf2dPoints_, ListOfSubdivParam_, orient_String:"x", Magnetiz_List:{0,0,0}, OptPar1_:0, OptPar2_:0, OptPar3_:0, OptPar4_:0]
:Arguments:     { N[xCoordin], N[lxWidth], N[ListOf2dPoints], N[ListOfSubdivParam], orient, N[Magnetiz[[1]]],N[Magnetiz[[2]]],N[Magnetiz[[3]]], OptPar1, OptPar2, OptPar3, OptPar4 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//void ArcMag P(( double,double,double, double,double, double,double, double, int, const char*, double,double,double ));
//
//:Begin:
//:Function:      ArcMag
//:Pattern:       radObjArcMag[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {phiminWidth_,phimaxWidth_}, HeightWidth_, SectN_, Orient_String:"z", Mv_List:{0,0,0}]
//:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[phiminWidth],N[phimaxWidth], N[HeightWidth], Round[SectN], Orient, N[Mv[[1]]],N[Mv[[2]]],N[Mv[[3]]] }
//:ArgumentTypes: { Real,Real,Real, Real,Real, Real,Real, Real, Integer, String, Real,Real,Real }
//:ReturnType:    Manual
//:End:


void ArcPolygon P(( ));

:Begin:
:Function:      ArcPolygon
:Pattern:       radObjArcPgnMag[{xCoordin_,yCoordin_}, Orient_String, ListOf2dPoints_, {phiminWidth_,phimaxWidth_}, SectN_, SymOrNoSym_String:"nosym", Magn_List:{0.,0.,0.}]
:Arguments:     { N[xCoordin],N[yCoordin], Orient, N[ListOf2dPoints], N[phiminWidth],N[phimaxWidth], Round[SectN], SymOrNoSym, N[Magn[[1]]],N[Magn[[2]]],N[Magn[[3]]] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//void CylMag P(( double,double,double, double, double, int, double,double,double, const char* ));

//:Begin:
//:Function:      CylMag
//:Pattern:       radObjCylMag[{xCoordin_,yCoordin_,zCoordin_}, Rad_, HeightWidth_, SectN_, Mv_List:{0,0,0}, Orient_String:"z"]
//:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[Rad], N[HeightWidth], Round[SectN], N[Mv[[1]]],N[Mv[[2]]],N[Mv[[3]]], Orient }
//:ArgumentTypes: { Real,Real,Real, Real, Real, Integer, Real,Real,Real, String }
//:ReturnType:    Manual
//:End:


void CylMag P(( double,double,double, double, double, int, const char*, double,double,double ));

:Begin:
:Function:      CylMag
:Pattern:       radObjCylMag[{xCoordin_,yCoordin_,zCoordin_}, Rad_, HeightWidth_, SectN_, Orient_String:"z", Mv_List:{0,0,0}]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[Rad], N[HeightWidth], Round[SectN], Orient, N[Mv[[1]]],N[Mv[[2]]],N[Mv[[3]]] }
:ArgumentTypes: { Real,Real,Real, Real, Real, Integer, String, Real,Real,Real }
:ReturnType:    Manual
:End:


void RecCur P(( double,double,double, double,double,double, double,double,double ));

:Begin:
:Function:      RecCur
:Pattern:       radObjRecCur[{xCoordin_,yCoordin_,zCoordin_}, {wxWidth_,wyWidth_,wzWidth_}, {jxCurrDens_,jyCurrDens_,jzCurrDens_}]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[wxWidth],N[wyWidth],N[wzWidth], N[jxCurrDens],N[jyCurrDens],N[jzCurrDens] }
:ArgumentTypes: { Real,Real,Real, Real,Real,Real, Real,Real,Real }
:ReturnType:    Manual
:End:


void ArcCur P(( double,double,double, double,double, double,double, double, int, double, const char*, const char* ));

:Begin:
:Function:      ArcCur
:Pattern:       radObjArcCur[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {phiminWidth_,phimaxWidth_}, HeightWidth_, SectN_, Jaz_, ManOrAuto_String:"man", Orient_String:"z"]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[phiminWidth],N[phimaxWidth], N[HeightWidth], Round[SectN], N[Jaz], ManOrAuto, Orient }
:ArgumentTypes: { Real,Real,Real, Real,Real, Real,Real, Real, Integer, Real, String, String }
:ReturnType:    Manual
:End:


void RaceTrack P(( double,double,double, double,double, double,double, double, int, double, const char*, const char* ));

:Begin:
:Function:      RaceTrack
:Pattern:       radObjRaceTrk[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {lxWidth_,lyWidth_}, HeightWidth_, SectN_, Jaz_, ManOrAuto_String:"man", Orient_String:"z"]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[lxWidth],N[lyWidth], N[HeightWidth], Round[SectN], N[Jaz], ManOrAuto, Orient }
:ArgumentTypes: { Real,Real,Real, Real,Real, Real,Real, Real, Integer, Real, String, String }
:ReturnType:    Manual
:End:


//:Begin:
//:Function:      RaceTrack
//:Pattern:       radObjRaceTrkCur[{xCoordin_,yCoordin_,zCoordin_}, {rminWidth_,rmaxWidth_}, {lxWidth_,lyWidth_}, HeightWidth_, SectN_, Jaz_, ManOrAuto_String:"man", Orient_String:"z"]
//:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[rminWidth],N[rmaxWidth], N[lxWidth],N[lyWidth], N[HeightWidth], Round[SectN], N[Jaz], ManOrAuto, Orient }
//:ArgumentTypes: { Real,Real,Real, Real,Real, Real,Real, Real, Integer, Real, String, String }
//:ReturnType:    Manual
//:End:


void FlmCur P(( ));

:Begin:
:Function:      FlmCur
:Pattern:       radObjFlmCur[ListOfPoints_, TotalCurrent_]
:Arguments:     { N[ListOfPoints], N[TotalCurrent] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void Rectngl P(( double,double,double, double,double ));

:Begin:
:Function:      Rectngl
:Pattern:       radFldFrcShpRtg[{xCoordin_,yCoordin_,zCoordin_},{wxWidth_,wyWidth_}]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[wxWidth],N[wyWidth] }
:ArgumentTypes: { Real,Real,Real, Real,Real }
:ReturnType:    Manual
:End:


void Group P(( int*, long ));

:Begin:
:Function:      Group
:Pattern:       radObjCnt[ArrayOfKeys_List]
:Arguments:     { ArrayOfKeys }
:ArgumentTypes: { IntegerList }
:ReturnType:    Manual
:End:


void AddToGroup P(( int, int*, long ));

:Begin:
:Function:      AddToGroup
:Pattern:       radObjAddToCnt[GroupKey_, ArrayOfKeys_List]
:Arguments:     { Round[GroupKey], ArrayOfKeys }
:ArgumentTypes: { Integer, IntegerList }
:ReturnType:    Manual
:End:


void OutGroupSubObjectKeys P(( int ));

:Begin:
:Function:      OutGroupSubObjectKeys
:Pattern:       radObjCntStuf[GroupKey_]
:Arguments:     { Round[GroupKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void BackgroundFieldSource P(( double,double,double ));

:Begin:
:Function:      BackgroundFieldSource
:Pattern:       radObjBckg[{bxField_,byField_,bzField_}]
:Arguments:     { N[bxField],N[byField],N[bzField] }
:ArgumentTypes: { Real,Real,Real }
:ReturnType:    Manual
:End:


void SubdivideElementG3D P(( ));

:Begin:
:Function:      SubdivideElementG3D
:Pattern:       radObjDivMag[ElemKey_, SubdivisionStructure_, OptAdditionalSpecifications_:0, OptPar1_:0, OptPar2_:0, OptPar3_:0]
:Arguments:     { Round[ElemKey], N[SubdivisionStructure], N[OptAdditionalSpecifications], OptPar1, OptPar2, OptPar3 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void CutElementG3D P(( ));

:Begin:
:Function:      CutElementG3D
:Pattern:       radObjCutMag[ElemKey_, PointOnThePlane_, NormalVector_, OptPar1_:0, OptPar2_:0]
:Arguments:     { Round[ElemKey], N[PointOnThePlane], N[NormalVector], OptPar1, OptPar2 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void SubdivideElementG3DByParPlanes P(( ));

:Begin:
:Function:      SubdivideElementG3DByParPlanes
:Pattern:       rAdObjDivByPlns[ElemKey_, SubdivisionStructure_, FrameFlag_String:"lab"]
:Arguments:     { Round[ElemKey], N[SubdivisionStructure], FrameFlag }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void DuplicateElementG3D P(( ));

:Begin:
:Function:      DuplicateElementG3D
:Pattern:       radObjDpl[ElemKey_, OptPar1_:0]
:Arguments:     { Round[ElemKey], OptPar1 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void CreateFromG3DObjectWithSymmetries P(( int ));

:Begin:
:Function:      CreateFromG3DObjectWithSymmetries
:Pattern:       rAdObjDplFrSym[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void NumberOfDegOfFreedom P(( int ));

:Begin:
:Function:      NumberOfDegOfFreedom
:Pattern:       radObjDegFre[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void MagnOfObj P(( int ));

:Begin:
:Function:      MagnOfObj
:Pattern:       radObjM[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void ObjField P(( int, const char* ));

:Begin:
:Function:      ObjField
:Pattern:       radObjCenFld[ElemKey_, FldCh_]
:Arguments:     { Round[ElemKey], FldCh }
:ArgumentTypes: { Integer, String }
:ReturnType:    Manual
:End:


void ScaleCurInObj P(( int,double ));

:Begin:
:Function:      ScaleCurInObj
:Pattern:       radObjScaleCur[ElemKey_, scaleCoef_]
:Arguments:     { Round[ElemKey], scaleCoef }
:ArgumentTypes: { Integer, Real }
:ReturnType:    Manual
:End:


void GeometricalVolume P(( int ));

:Begin:
:Function:      GeometricalVolume
:Pattern:       rAdObjGeoVol[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:

:Begin:
:Function:      GeometricalVolume
:Pattern:       radObjGeoVol[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void GeometricalLimits P(( int ));

:Begin:
:Function:      GeometricalLimits
:Pattern:       radObjGeoLim[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void FldCmpMetForSubdRecMag P(( int, int, int ));

:Begin:
:Function:      FldCmpMetForSubdRecMag
:Pattern:       rAdFldCmpMetNxNyNz[ElemKey_, MethNo_, SubLevel_:0]
:Arguments:     { Round[ElemKey], Round[MethNo], Round[SubLevel] }
:ArgumentTypes: { Integer, Integer, Integer }
:ReturnType:    Manual
:End:


void SetLocMgnInSbdRecMag P(( ));

:Begin:
:Function:      SetLocMgnInSbdRecMag
:Pattern:       rAdObjSetLocMgn[ElemKey_, ListOfMagnetiz_]
:Arguments:     { Round[ElemKey], ListOfMagnetiz }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//-------------------------------------------------------------------------


void Translation P(( double,double,double ));

:Begin:
:Function:      Translation
:Pattern:       radTrfTrsl[{vxVector_,vyVector_,vzVector_}]
:Arguments:     { N[vxVector],N[vyVector],N[vzVector] }
:ArgumentTypes: { Real,Real,Real }
:ReturnType:    Manual
:End:


void Rotation P(( double,double,double, double,double,double, double ));

:Begin:
:Function:      Rotation
:Pattern:       radTrfRot[{xCoordin_,yCoordin_,zCoordin_}, {vxVector_,vyVector_,vzVector_}, Angle_]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[vxVector],N[vyVector],N[vzVector], N[Angle] }
:ArgumentTypes: { Real,Real,Real, Real,Real,Real, Real }
:ReturnType:    Manual
:End:


void PlaneSym P(( double,double,double, double,double,double ));

:Begin:
:Function:      PlaneSym
:Pattern:       radTrfPlSym[{xCoordin_,yCoordin_,zCoordin_}, {nxVector_,nyVector_,nzVector_}]
:Arguments:     { N[xCoordin],N[yCoordin],N[zCoordin], N[nxVector],N[nyVector],N[nzVector] }
:ArgumentTypes: { Real,Real,Real, Real,Real,Real }
:ReturnType:    Manual
:End:


void FieldInversion P(( ));

:Begin:
:Function:      FieldInversion
:Pattern:       radTrfInv[]
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Manual
:End:


void TransformObject P(( int, int ));

:Begin:
:Function:      TransformObject
:Pattern:       radTrfOrnt[g3dElemKey_, TransElemKey_]
:Arguments:     { Round[g3dElemKey], Round[TransElemKey] }
:ArgumentTypes: { Integer, Integer }
:ReturnType:    Manual
:End:


void ApplySymmetry P(( int, int, int ));

:Begin:
:Function:      ApplySymmetry
:Pattern:       radTrfMlt[g3dElemKey_, TransElemKey_, Multiplicity_]
:Arguments:     { Round[g3dElemKey], Round[TransElemKey], Round[Multiplicity] }
:ArgumentTypes: { Integer, Integer, Integer }
:ReturnType:    Manual
:End:


void CombineTransformLeft P(( int, int ));

:Begin:
:Function:      CombineTransformLeft
:Pattern:       radTrfCmbL[ThisElemKey_, AnotherElemKey_]
:Arguments:     { Round[ThisElemKey], Round[AnotherElemKey] }
:ArgumentTypes: { Integer, Integer }
:ReturnType:    Manual
:End:


void CombineTransformRight P(( int, int ));

:Begin:
:Function:      CombineTransformRight
:Pattern:       radTrfCmbR[ThisElemKey_, AnotherElemKey_]
:Arguments:     { Round[ThisElemKey], Round[AnotherElemKey] }
:ArgumentTypes: { Integer, Integer }
:ReturnType:    Manual
:End:


//-------------------------------------------------------------------------


void LinearMaterial P(( double,double, double,double,double ));

:Begin:
:Function:      LinearMaterial
:Pattern:       radMatLin[{ksiparMagneticStuff_,ksiperMagneticStuff_}, {mrxMagneticStuff_,mryMagneticStuff_,mrzMagneticStuff_}]
:Arguments:     { N[ksiparMagneticStuff],N[ksiperMagneticStuff], N[mrxMagneticStuff],N[mryMagneticStuff],N[mrzMagneticStuff] }
:ArgumentTypes: { Real,Real, Real,Real,Real }
:ReturnType:    Manual
:End:


void LinearMaterial2 P(( double,double, double ));

:Begin:
:Function:      LinearMaterial2
:Pattern:       radMatLin[{ksiparMagneticStuff_,ksiperMagneticStuff_}, RemanentMagnetizMagnitude_]
:Arguments:     { N[ksiparMagneticStuff],N[ksiperMagneticStuff], N[RemanentMagnetizMagnitude] }
:ArgumentTypes: { Real,Real, Real }
:ReturnType:    Manual
:End:


void NonlinearIsotropMaterial P(( double,double,double, double,double,double ));

:Begin:
:Function:      NonlinearIsotropMaterial
:Pattern:       radMatSatIso[{ms1MagneticStuff_,ms2MagneticStuff_,ms3MagneticStuff_},{ks1MagneticStuff_,ks2MagneticStuff_,ks3MagneticStuff_}]
:Arguments:     { N[ms1MagneticStuff],N[ms2MagneticStuff],N[ms3MagneticStuff], N[ks1MagneticStuff],N[ks2MagneticStuff],N[ks3MagneticStuff] }
:ArgumentTypes: { Real,Real,Real, Real,Real,Real }
:ReturnType:    Manual
:End:


void NonlinearIsotropMaterial2 P(( double,double, double,double, double,double ));

:Begin:
:Function:      NonlinearIsotropMaterial2
:Pattern:       radMatSatIso[{ks1MagneticStuff_, ms1MagneticStuff_}, {ks2MagneticStuff_, ms2MagneticStuff_}, {ks3MagneticStuff_, ms3MagneticStuff_}]
:Arguments:     { N[ks1MagneticStuff],N[ms1MagneticStuff], N[ks2MagneticStuff],N[ms2MagneticStuff], N[ks3MagneticStuff],N[ms3MagneticStuff] }
:ArgumentTypes: { Real,Real, Real,Real, Real,Real }
:ReturnType:    Manual
:End:


void NonlinearIsotropMaterial3 P(( ));

:Begin:
:Function:      NonlinearIsotropMaterial3
:Pattern:       radMatSatIso[ListOfPairsHM_]
:Arguments:     { N[ListOfPairsHM] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void NonlinearLaminatedMaterialML P(( ));

:Begin:
:Function:      NonlinearLaminatedMaterialML
:Pattern:       radMatSatLam[MagData_, LamCoef_, LamNormal_:0]
:Arguments:     { N[MagData], N[LamCoef], N[LamNormal] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void NonlinearAnisotropMaterial P(( ));

:Begin:
:Function:      NonlinearAnisotropMaterial
:Pattern:       radMatSatAniso[MagDataPara_, MagDataPerp_]
:Arguments:     { N[MagDataPara], N[MagDataPerp] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void ApplyMaterial P(( int, int ));

:Begin:
:Function:      ApplyMaterial
:Pattern:       radMatApl[g3dRelaxElemKey_, MaterElemKey_]
:Arguments:     { Round[g3dRelaxElemKey], Round[MaterElemKey] }
:ArgumentTypes: { Integer, Integer }
:ReturnType:    Manual
:End:


void MvsH P(( int, const char*, double,double,double ));

:Begin:
:Function:      MvsH
:Pattern:       radMatMvsH[g3dRelaxOrMaterElemKey_, MagnChar_String, {hxMagneticStuff_,hyMagneticStuff_,hzMagneticStuff_}]
:Arguments:     { Round[g3dRelaxOrMaterElemKey], MagnChar, N[hxMagneticStuff],N[hyMagneticStuff],N[hzMagneticStuff] }
:ArgumentTypes: { Integer, String, Real,Real,Real }
:ReturnType:    Manual
:End:


//-------------------------------------------------------------------------


void PreRelax P(( int, int ));

:Begin:
:Function:      PreRelax
:Pattern:       radRlxPre[ElemKey_, SrcElemKey_:0]
:Arguments:     { Round[ElemKey], Round[SrcElemKey] }
:ArgumentTypes: { Integer, Integer }
:ReturnType:    Manual
:End:


void ShowInteractMatrix P(( int ));

:Begin:
:Function:      ShowInteractMatrix
:Pattern:       rAdRlxOutIntrcMatr[InteractElemKey_]
:Arguments:     { Round[InteractElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void ShowInteractVector P(( int, const char* ));

:Begin:
:Function:      ShowInteractVector
:Pattern:       rAdRlxOutIntrcVect[InteractElemKey_, InteractFieldID_String]
:Arguments:     { Round[InteractElemKey], InteractFieldID }
:ArgumentTypes: { Integer, String }
:ReturnType:    Manual
:End:


void ManualRelax P(( int, int, int, double ));

:Begin:
:Function:      ManualRelax
:Pattern:       radRlxMan[InteractElemKey_, MethNo_, IterNum_, RelaxPar_]
:Arguments:     { Round[InteractElemKey], Round[MethNo], Round[IterNum], N[RelaxPar] }
:ArgumentTypes: { Integer, Integer, Integer, Real }
:ReturnType:    Manual
:End:


//void AutoRelax P(( int, double, int, int ));
void AutoRelax P(( ));

:Begin:
:Function:      AutoRelax
:Pattern:       radRlxAuto[InteractElemKey_, Prec_, MaxIterNum_, MethNo_:-1, OptPar1_:-1]
:Arguments:     { Round[InteractElemKey], N[Prec], Round[MaxIterNum], MethNo, OptPar1 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void UpdateSourcesForRelax P(( int ));

:Begin:
:Function:      UpdateSourcesForRelax
:Pattern:       radRlxUpdSrc[InteractElemKey_]
:Arguments:     { Round[InteractElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void SolveGen P(( int, double, int, int ));

:Begin:
:Function:      SolveGen
:Pattern:       radSolve[ObjKey_, Prec_, MaxIterNum_, MethNo_:4]
:Arguments:     { Round[ObjKey], N[Prec], Round[MaxIterNum], Round[MethNo] }
:ArgumentTypes: { Integer, Real, Integer, Integer }
:ReturnType:    Manual
:End:


void CompCriterium P(( double, double, double, double, double,double ));

:Begin:
:Function:      CompCriterium
:Pattern:       radFldCmpCrt[AbsPrecB_, AbsPrecA_, AbsPrecBInt_, AbsPrecFrc_:1., AbsPrecTrjCrd_:-1, AbsPrecTrjAng_:-1]
:Arguments:     { N[AbsPrecB], N[AbsPrecA], N[AbsPrecBInt], N[AbsPrecFrc], N[AbsPrecTrjCrd],N[AbsPrecTrjAng] }
:ArgumentTypes: { Real, Real, Real, Real, Real,Real }
:ReturnType:    Manual
:End:


void CompPrecision P(( ));

:Begin:
:Function:      CompPrecision
:Pattern:       radFldCmpPrc[Prc1_, Prc2_:0, Prc3_:0, Prc4_:0, Prc5_:0, Prc6_:0, Prc7_:0, Prc8_:0, Prc9_:0]
:Arguments:     { Prc1, Prc2, Prc3, Prc4, Prc5, Prc6, Prc7, Prc8, Prc9 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void MultipoleThresholds P(( double, double, double, double )); // Maybe to be removed later

:Begin:
:Function:      MultipoleThresholds
:Pattern:       rAdFldMltpolThr[{a0Threshold_,a1Threshold_,a2Threshold_,a3Threshold_}]
:Arguments:     { N[a0Threshold],N[a1Threshold],N[a2Threshold],N[a3Threshold] }
:ArgumentTypes: { Real,Real,Real,Real }
:ReturnType:    Manual
:End:


void Field P(( int, const char*, double,double,double, double,double,double, int, const char*, double ));

//:Begin:
//:Function:      Field
//:Pattern:       radFld[ElemKey_, FieldChar_String, {xCoordin_,yCoordin_,zCoordin_}, x2Coordin_:10.^23,y2Coordin_:10.^23,z2Coordin_:10.^23, np_:1, ShowArgFlag_String:"noarg", strtarg_:0.]
//:Arguments:     { Round[ElemKey], FieldChar, N[xCoordin],N[yCoordin],N[zCoordin], N[x2Coordin],N[y2Coordin],N[z2Coordin], Round[np], ShowArgFlag, N[strtarg]}
//:ArgumentTypes: { Integer, String, Real,Real,Real, Real,Real,Real, Integer, String, Real }
//:ReturnType:    Manual
//:End:


:Begin:
:Function:      Field
:Pattern:       radFldLst[ElemKey_, FieldChar_String, {xCoordin_,yCoordin_,zCoordin_}, {x2Coordin_,y2Coordin_,z2Coordin_}, np_:1, ShowArgFlag_String:"noarg", strtarg_:0.]
:Arguments:     { Round[ElemKey], FieldChar, N[xCoordin],N[yCoordin],N[zCoordin], N[x2Coordin],N[y2Coordin],N[z2Coordin], Round[np], ShowArgFlag, N[strtarg]}
:ArgumentTypes: { Integer, String, Real,Real,Real, Real,Real,Real, Integer, String, Real }
:ReturnType:    Manual
:End:


void ShimSignature P(( int, const char*, double,double,double, double,double,double, double,double,double, int, double,double,double ));

:Begin:
:Function:      ShimSignature
:Pattern:       radFldShimSig[ElemKey_, FieldChar_String, {dxCoordin_,dyCoordin_,dzCoordin_}, {xCoordin_,yCoordin_,zCoordin_}, {x2Coordin_,y2Coordin_,z2Coordin_}, np_, {vxCoordin_:0.,vyCoordin_:0.,vzCoordin_:0.}]
:Arguments:     { Round[ElemKey], FieldChar, N[dxCoordin],N[dyCoordin],N[dzCoordin], N[xCoordin],N[yCoordin],N[zCoordin], N[x2Coordin],N[y2Coordin],N[z2Coordin], Round[np], N[vxCoordin],N[vyCoordin],N[vzCoordin]}
:ArgumentTypes: { Integer, String, Real,Real,Real, Real,Real,Real, Real,Real,Real, Integer, Real,Real,Real }
:ReturnType:    Manual
:End:


void FieldArbitraryPointsStruct P(( ));

//:Begin:
//:Function:      FieldArbitraryPointsStruct
//:Pattern:       radFld[ElemKey_, FieldChar_String, PointsStructure_]
//:Arguments:     { Round[ElemKey], FieldChar, N[PointsStructure] }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:

:Begin:
:Function:      FieldArbitraryPointsStruct
:Pattern:       radFld[ElemKey_, FieldChar_String, PointsStructure_]
:Arguments:     { Round[ElemKey], FieldChar, N[PointsStructure] }
:ArgumentTypes: { Integer, String, Manual }
:ReturnType:    Manual
:End:


void FieldInt P(( int, const char*, const char*, double,double,double, double,double,double ));

:Begin:
:Function:      FieldInt
:Pattern:       radFldInt[ElemKey_, CondChar_String, FieldIntChar_String, {x1Coordin_,y1Coordin_,z1Coordin_}, {x2Coordin_,y2Coordin_,z2Coordin_}]
:Arguments:     { Round[ElemKey], CondChar, FieldIntChar, N[x1Coordin],N[y1Coordin],N[z1Coordin], N[x2Coordin],N[y2Coordin],N[z2Coordin] }
:ArgumentTypes: { Integer, String, String, Real,Real,Real, Real,Real,Real }
:ReturnType:    Manual
:End:


void FieldForce P(( int, int ));

:Begin:
:Function:      FieldForce
:Pattern:       radFldFrc[ObjElemKey_, ShapeElemKey_]
:Arguments:     { Round[ObjElemKey], Round[ShapeElemKey] }
:ArgumentTypes: { Integer, Integer }
:ReturnType:    Manual
:End:


void FieldEnergy P(( int, int, int,int,int ));

:Begin:
:Function:      FieldEnergy
:Pattern:       radFldEnr[DestObjElemKey_, SourceObjElemKey_, SubdivParam_List:{0,0,0}]
:Arguments:     { Round[DestObjElemKey], Round[SourceObjElemKey], Round[SubdivParam[[1]]],Round[SubdivParam[[2]]],Round[SubdivParam[[3]]] }
:ArgumentTypes: { Integer, Integer, Integer,Integer,Integer }
:ReturnType:    Manual
:End:


void FieldForceThroughEnergy P(( int, int, const char*, int,int,int ));

:Begin:
:Function:      FieldForceThroughEnergy
:Pattern:       radFldEnrFrc[DestObjElemKey_, SourceObjElemKey_, ComponIDChar_String, SubdivParam_List:{0,0,0}]
:Arguments:     { Round[DestObjElemKey], Round[SourceObjElemKey], ComponIDChar, Round[SubdivParam[[1]]],Round[SubdivParam[[2]]],Round[SubdivParam[[3]]] }
:ArgumentTypes: { Integer, Integer, String, Integer,Integer,Integer }
:ReturnType:    Manual
:End:


void FieldTorqueThroughEnergy P(( int, int, const char*, double,double,double, int,int,int ));

:Begin:
:Function:      FieldTorqueThroughEnergy
:Pattern:       radFldEnrTrq[DestObjElemKey_, SourceObjElemKey_, ComponIDChar_String, TorqueCentrPo_List:{0.,0.,0.}, SubdivParam_List:{0,0,0}]
:Arguments:     { Round[DestObjElemKey], Round[SourceObjElemKey], ComponIDChar, N[TorqueCentrPo[[1]]],N[TorqueCentrPo[[2]]],N[TorqueCentrPo[[3]]], Round[SubdivParam[[1]]],Round[SubdivParam[[2]]],Round[SubdivParam[[3]]] }
:ArgumentTypes: { Integer, Integer, String, Real,Real,Real, Integer,Integer,Integer }
:ReturnType:    Manual
:End:


void ParticleTrajectory P(( int, double, double,double,double,double, double,double, int ));

:Begin:
:Function:      ParticleTrajectory
:Pattern:       radFldPtcTrj[ElemKey_,Energy_,{x0Coordin_,dxdy0_,z0Coordin_,dzdy0_},{y0Coordin_,y1Coordin_},np_]
:Arguments:     { Round[ElemKey], N[Energy], N[x0Coordin],N[dxdy0],N[z0Coordin],N[dzdy0], N[y0Coordin],N[y1Coordin], Round[np] }
:ArgumentTypes: { Integer, Real, Real,Real,Real,Real, Real,Real, Integer }
:ReturnType:    Manual
:End:


void FocusingPotential P(( int, double,double,double, double,double,double, int ));

:Begin:
:Function:      FocusingPotential
:Pattern:       radFldFocPot[ElemKey_,{xxStart_,yyStart_,zzStart_},{xxFin_,yyFin_,zzFin_},NumPo_]
:Arguments:     { Round[ElemKey], N[xxStart],N[yyStart],N[zzStart], N[xxFin],N[yyFin],N[zzFin], Round[NumPo] }
:ArgumentTypes: { Integer, Real,Real,Real, Real,Real,Real, Integer }
:ReturnType:    Manual
:End:


//void FocusingKickPer P(( int, double,double,double, double,double,double, double,int, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double ));
void FocusingKickPer P(( int, double,double,double, double,double,double, double,double, double,double,double, double,int,double,int, const char*, int,int,double,double, const char*, double, const char* ));

:Begin:
:Function:      FocusingKickPer
:Pattern:       radFldFocKickPer[ElemKey_,{x0_,y0_,z0_},{nsx_,nsy_,nsz_},per_,nper_,{n1x_,n1y_,n1z_},r1_,np1_,r2_,np2_,comment_String:"",{nh_:1,ns_:8,d1_:0,d2_:0},kickUnit_String:"T2m2",en_:1,format_String:"fix"]
:Arguments:     { Round[ElemKey], N[x0],N[y0],N[z0], N[nsx],N[nsy],N[nsz], N[per],N[nper], N[n1x],N[n1y],N[n1z], N[r1],Round[np1],N[r2],Round[np2], comment, Round[nh],Round[ns],N[d1],N[d2], kickUnit, N[en], format }
:ArgumentTypes: { Integer, Real,Real,Real, Real,Real,Real, Real,Real, Real,Real,Real, Real,Integer,Real,Integer, String, Integer,Integer,Real,Real, String, Real, String }
:ReturnType:    Manual
:End:


//void FocusingKick P(( int, double,double,double, double,double,double, double*,long,int, double,double,double, double,int,double,int, const char*, double,double ));
void FocusingKickML P(( ));

:Begin:
:Function:      FocusingKickML
:Pattern:       radFldFocKick[ElemKey_,P1_List,Nlong_List,ArrayOfLongPos_List,nps_,Ntr1_List,r1_,np1_,r2_,np2_,comment_String:"",{d1_:0,d2_:0}]
:Arguments:     { Round[ElemKey], N[P1], N[Nlong], N[ArrayOfLongPos],Round[nps], N[Ntr1], N[r1],Round[np1],N[r2],Round[np2], comment, N[d1],N[d2] }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void PhysicalUnits P(( ));

:Begin:
:Function:      PhysicalUnits
:Pattern:       radFldUnits[]
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Manual
:End:


//void OffsetForConvergence P(( char*, double,double,double, double,double,double ));

//:Begin:
//:Function:      OffsetForConvergence
//:Pattern:       radFldOfst[PointsOrDims_String, {axReg_,ayReg_,azReg_}, {axRnd_,ayRnd_,azRnd_}]
//:Arguments:     { PointsOrDims, N[axReg],N[ayReg],N[azReg], N[axRnd],N[ayRnd],N[azRnd] }
//:ArgumentTypes: { String, Real,Real,Real, Real,Real,Real }
//:ReturnType:    Manual
//:End:


void TolForConvergence P(( double, double, double ));

:Begin:
:Function:      TolForConvergence
:Pattern:       radFldLenTol[AbsRandMagnitude_,RelRandMagnitude_,ZeroRandMagnitude_:0]
:Arguments:     { N[AbsRandMagnitude], N[RelRandMagnitude], N[ZeroRandMagnitude] }
:ArgumentTypes: { Real, Real, Real }
:ReturnType:    Manual
:End:


void RandomizationOnOrOff P(( const char* ));

:Begin:
:Function:      RandomizationOnOrOff
:Pattern:       radFldLenRndSw[OnOrOff_String]
:Arguments:     { OnOrOff }
:ArgumentTypes: { String }
:ReturnType:    Manual
:End:


//-------------------------------------------------------------------------


void ApplyDrawAttrToElem P(( int, double,double,double, double ));

:Begin:
:Function:      ApplyDrawAttrToElem
:Pattern:       radObjDrwAtr[ElemKey_, {rColorAttrib_,gColorAttrib_,bColorAttrib_}, LineThickness_: -1.]
:Arguments:     { Round[ElemKey], N[rColorAttrib],N[gColorAttrib],N[bColorAttrib], N[LineThickness] }
:ArgumentTypes: { Integer, Real,Real,Real, Real }
:ReturnType:    Manual
:End:


//void ApplyColorToElem P(( int, double,double,double ));

//:Begin:
//:Function:      ApplyColorToElem
//:Pattern:       radObjCol[ElemKey_, {r_,g_,b_}]
//:Arguments:     { Round[ElemKey], N[r],N[g],N[b] }
//:ArgumentTypes: { Integer, Real,Real,Real }
//:ReturnType:    Manual
//:End:


void RemoveDrawAttrFromElem P(( int ));

:Begin:
:Function:      RemoveDrawAttrFromElem
:Pattern:       rAdObjDelDrwAtr[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


//void GraphicsForElemWithSymChilds P(( int ));

//:Begin:
//:Function:      GraphicsForElemWithSymChilds
//:Pattern:       radObjDrw[ElemKey_]
//:Arguments:     { Round[ElemKey] }
//:ArgumentTypes: { Integer }
//:ReturnType:    Manual
//:End:


void GraphicsForElemWithSymChildsExt P(( ));

:Begin:
:Function:      GraphicsForElemWithSymChildsExt
:Pattern:       radObjDrw[ElemKey_, OptPar1_:0]
:Arguments:     { Round[ElemKey], OptPar1 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//void GraphicsForAllWithSymChilds P(( ));

//:Begin:
//:Function:      GraphicsForAllWithSymChilds
//:Pattern:       radObjDrwAll
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


void GraphicsForElemWithoutSymChilds P(( int ));

:Begin:
:Function:      GraphicsForElemWithoutSymChilds
:Pattern:       rAdObjDrwWithoutTrfMlt[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


//void GraphicsForAllWithoutSymChilds P(( ));

//:Begin:
//:Function:      GraphicsForAllWithoutSymChilds
//:Pattern:       radObjDrwAllWithoutTrfMlt
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


//void QuickDraw3D_Viewer P(( ));

//:Begin:
//:Function:      QuickDraw3D_Viewer
//:Pattern:       radObjDrwQD3D[ElemKey_, OptPar1_:0, OptPar2_:0, OptPar3_:0]
//:Arguments:     { Round[ElemKey], OptPar1, OptPar2, OptPar3 }
//:ArgumentTypes: { Manual }
//:ReturnType:    Manual
//:End:


void OpenGL_3D_Viewer P(( ));

:Begin:
:Function:      OpenGL_3D_Viewer
:Pattern:       radObjDrwOpenGL[ElemKey_, OptPar1_:0, OptPar2_:0, OptPar3_:0]
:Arguments:     { Round[ElemKey], OptPar1, OptPar2, OptPar3 }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


//-------------------------------------------------------------------------


void DeleteElement P(( int ));

:Begin:
:Function:      DeleteElement
:Pattern:       radUtiDel[ElemKey_]
:Arguments:     { Round[ElemKey] }
:ArgumentTypes: { Integer }
:ReturnType:    Manual
:End:


void DeleteAllElements1 P(( ));

:Begin:
:Function:      DeleteAllElements1
:Pattern:       radUtiDelAll[]
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Manual
:End:


void DeleteAllElements2 P(( ));

:Begin:
:Function:      DeleteAllElements2
:Pattern:       rAdUtiDelAll[]
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Manual
:End:


//void DumpElem P(( int ));

//:Begin:
//:Function:      DumpElem
//:Pattern:       radUtiDmp[ElemKey_, OutFormat_String:"asc"]
//:Arguments:     { Round[ElemKey], OutFormat }
//:ArgumentTypes: { Integer, String }
//:ReturnType:    Manual
//:End:

void DumpElem P(( ));

:Begin:
:Function:      DumpElem
:Pattern:       radUtiDmp[ElemKey_, OutFormat_String:"asc"]
:Arguments:     { Round[ElemKey], OutFormat }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:

//void GenDump P(( ));

//:Begin:
//:Function:      GenDump
//:Pattern:       radUtiDmpAll
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


void DumpElemParse P(( ));

:Begin:
:Function:      DumpElemParse
:Pattern:       radUtiDmpPrs[bstr_]
:Arguments:     { bstr }
:ArgumentTypes: { Manual }
:ReturnType:    Manual
:End:


void RadiaVersion P(());

:Begin:
:Function:      RadiaVersion
:Pattern:       radUtiVer[]
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Manual
:End:


//void OutCommandsInfo P(( ));

//:Begin:
//:Function:      OutCommandsInfo
//:Pattern:       radUtiInfo
//:Arguments:     { }
//:ArgumentTypes: { }
//:ReturnType:    Manual
//:End:


void ReturnInput P(( double, int ));

:Begin:
:Function:      ReturnInput
:Pattern:       rAdUtiRetInp[Input_, NumTimes_]
:Arguments:     { N[Input], Round[NumTimes] }
:ArgumentTypes: { Real, Integer }
:ReturnType:    Manual
:End:


void MemAllocMethForIntrctMatr P(( const char* ));

:Begin:
:Function:      MemAllocMethForIntrctMatr
:Pattern:       rAdUtiRlxMemAllocMet[TotOrParts_String]
:Arguments:     { TotOrParts }
:ArgumentTypes: { String }
:ReturnType:    Manual
:End:


//------------ P ELLEAUME -------------------------------------------------


void StartProf P(( int, int, int ));

:Begin:
:Function:      StartProf
:Pattern:       rAdUtiStrtProf[Flag_, NumFunc_, Depth_]
:Arguments:     { Round[Flag], Round[NumFunc], Round[Depth] }
:ArgumentTypes: { Integer, Integer, Integer }
:ReturnType:    Manual
:End:


void StopProf P(());

:Begin:
:Function:      StopProf
:Pattern:       rAdUtiStopProf
:Arguments:     { }
:ArgumentTypes: { }
:ReturnType:    Manual
:End:


void InterruptTime P(( double ));

:Begin:
:Function:       InterruptTime
:Pattern:        radUtiIntrptTim[tTimeQuanta_:1.0]
:Arguments:      { N[tTimeQuanta] }
:ArgumentTypes:  { Real }
:ReturnType:     Manual
:End:


//-------------------------------------------------------------------------


//:Evaluate:      End[]


//-------------------------------------------------------------------------


//:Evaluate:      EndPackage[]

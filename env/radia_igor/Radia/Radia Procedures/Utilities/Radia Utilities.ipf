#pragma rtGlobals=1		// Use modern global access method

//+++++++++++++++++++++++++++++++++++++++
//
//General Utilities
//
//+++++++++++++++++++++++++++++++++++++++
//Show Help Topic from non-formatted Help File notebook.
//+++++++++++++++++++++++++++++++++++++++
function RadUtiShowHelpTopic(TopicStr)
string TopicStr

string HelpNotebookName = "Radia Help.ifn"
opennotebook/R/P=RadUtiPath/N=RadiaHelpNotebook HelpNotebookName
notebook RadiaHelpNotebook selection={startOfFile, startOfFile}
notebook RadiaHelpNotebook findText={TopicStr, 1}
MoveWindow 150, 150, 650, 550
end

//+++++++++++++++++++++++++++++++++++++++
//Show Radia Example procedure file.
//+++++++++++++++++++++++++++++++++++++++
function RadUtiShowProcExample(fname)
string fname
pathinfo RadUtiPath
newpath/O/Q RadUtiPathEx, S_path+"Radia Procedures:Examples"
opennotebook/Z/P=RadUtiPathEx/N=RadiaExample/T="IPRC" fname
//open/Z/R/P=RadUtiPathEx/N=RadiaExample/T="IPRC" fname
end

//+++++++++++++++++++++++++++++++++++++++
//
//+++++++++++++++++++++++++++++++++++++++
proc RadUtiKillAllGraphs()
DoWindow/K Graph0
DoWindow/K Graph1
DoWindow/K Graph2
DoWindow/K Graph3
DoWindow/K Graph4
DoWindow/K Graph5
DoWindow/K Graph6
DoWindow/K Graph7
DoWindow/K Graph8
DoWindow/K Graph9
DoWindow/K Graph10
DoWindow/K Graph11
DoWindow/K Graph12
DoWindow/K Graph13
DoWindow/K Graph14
DoWindow/K Graph15
DoWindow/K Graph16
DoWindow/K Graph17
DoWindow/K Graph18
DoWindow/K Graph19
DoWindow/K Graph20
DoWindow/K Graph21
DoWindow/K Graph22
DoWindow/K Graph23
DoWindow/K Graph24
DoWindow/K Graph25
DoWindow/K Graph26
DoWindow/K Graph27
DoWindow/K Graph28
DoWindow/K Graph29
DoWindow/K Graph30
DoWindow/K Graph31
DoWindow/K Graph32
DoWindow/K Graph33
DoWindow/K Graph34
DoWindow/K Graph35
DoWindow/K Graph36
DoWindow/K Graph37
DoWindow/K Graph38
DoWindow/K Graph39
DoWindow/K Graph40
DoWindow/K Graph41
DoWindow/K Graph42
DoWindow/K Graph43
DoWindow/K Graph44
DoWindow/K Graph45
DoWindow/K Graph46
DoWindow/K Graph47
DoWindow/K Graph48
DoWindow/K Graph49
DoWindow/K Graph50
// Make more smart
end

//+++++++++++++++++++++++++++++++++++++++
//
//+++++++++++++++++++++++++++++++++++++++
proc RadUtiKillExamWindows(DelayTime)
Variable DelayTime = 3 // Time delay before killing example graphs
DoUpdate; DoWindow/K srwHelp; Sleep/S DelayTime; RadUtiKillAllGraphs()
end

//+++++++++++++++++++++++++++++++++++++++
//Runs all Examples one-by-one
//+++++++++++++++++++++++++++++++++++++++
proc RadUtiAllExam()

RadiaExample1(); RadUtiKillExamWindows(2)
RadiaExample2(); RadUtiKillExamWindows(2)
RadiaExample3(); RadUtiKillExamWindows(2)
RadiaExample4(); RadUtiKillExamWindows(2)

//add next example here

end

//+++++++++++++++++++++++++++++++++++++++
//Auxiliary function to facilitate calculation of field integrals
//+++++++++++++++++++++++++++++++++++++++
function radFldIntCmp(Obj, CmpIb, CmpDir, xc, yc, zc)
variable Obj
string CmpIb, CmpDir
variable xc, yc, zc

make/O DummyIB
make/O FldIntCmpP1={xc, yc, zc}, FldIntCmpP2={xc, yc, zc}
variable ic = 0
if((cmpstr(CmpDir, "x") == 0) %| (cmpstr(CmpDir, "X") == 0))
	ic = 0
else
	if((cmpstr(CmpDir, "y") == 0) %| (cmpstr(CmpDir, "Y") == 0))
		ic = 1
	else
		if((cmpstr(CmpDir, "z") == 0) %| (cmpstr(CmpDir, "Z") == 0))
			ic = 2
		else
			abort "The direction string should be \"x\", \"y\", \"z\"."
		endif
	endif
endif
FldIntCmpP1[ic] -= 1000; FldIntCmpP2[ic] += 1000
variable ResIb = radFldInt(DummyIB, Obj, "inf", CmpIb, FldIntCmpP1,FldIntCmpP2)
killwaves/Z DummyIB, FldIntCmpP1, FldIntCmpP2
return ResIb
end


//+++++++++++++++++++++++++++++++++++++++
//
// Menu Radia
//
//+++++++++++++++++++++++++++++++++++++++
Menu "Radia"

"Help", RadUtiShowHelpTopic("About Radia     ")
"Initialize", RadInit(1)
"-"
"Example # 1: Uniformly Magnetized Cube", RadiaExample1()
"Example # 2: Superconducting Wiggler [Coils]", RadiaExample2()
"Example # 3: Generic Hybrid Undulator", RadiaExample3()
"Example # 4: Uniformly Magnetized Polyhedron", RadiaExample4()
"-"
"Run All Radia Examples",RadUtiAllExam()

end  // Menu Radia


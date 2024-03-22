## declaration of global functions and methods...

DeclareGlobalFunction("cos2piover");
DeclareGlobalFunction("sin2piover");
DeclareGlobalFunction("RotationMatrixDim2");
DeclareGlobalFunction("RotationMatrixDim3");

DeclareGlobalFunction("DimFixed");
DeclareGlobalFunction("DimFixedByGroup");
DeclareGlobalFunction("ConfigurationSpaceChar");

DeclareAttribute("SRC_QQ",IsLagSymmetryGroup);
DeclareGlobalFunction("DimChiKerT");

DeclareGlobalFunction("IsotropySubgroups");
DeclareGlobalFunction("MaximalIsotropies");
DeclareGlobalFunction("NumberOfReflections");

DeclareGlobalFunction("IsTypeRDirection");
DeclareGlobalFunction("TypeRDirections");

DeclareGlobalFunction("TransitiveDecomposition");


## methods:
IsValidLSG:=NewProperty("IsValidLSG",IsLagSymmetryGroup);

##
HasAlwaysCollisions:=NewProperty("HasAlwaysCollisions",IsLagSymmetryGroup);
IsCoercive:=NewProperty("IsCoercive",IsLagSymmetryGroup);


##########################################################################
IsCoerciveCHARS:=NewProperty("IsCoerciveCHARS",IsLagSymmetryGroupCHARS);
HasAlwaysCollisionsCHARS:=NewProperty("HasAlwaysCollisionsCHARS",IsLagSymmetryGroupCHARS);
HasAlwaysCollisionsCHARStemp:=NewProperty("HasAlwaysCollisionsCHARStemp",IsLagSymmetryGroupCHARS);
DeclareGlobalFunction("DimChiKerTCHARS");
##########################################################################




IsTypeR:=NewProperty("IsTypeR",IsLagSymmetryGroup); 

IsTransitiveLSG:=NewProperty("IsTransitiveLSG",IsLagSymmetryGroup); 

##
IsRedundant:=NewProperty("IsRedundant",IsLagSymmetryGroup);
DeclareAttribute("ActionType",IsLagSymmetryGroup);
DeclareAttribute("GroupOrder",IsLagSymmetryGroup);
DeclareAttribute("KernelTauOrder",IsLagSymmetryGroup);


